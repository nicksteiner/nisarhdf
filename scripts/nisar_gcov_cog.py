#!/usr/bin/env python3
"""
Convert NISAR L2 GCOV granules into Cloud-Optimized GeoTIFFs and PySTAC items.
Georeferencing and metadata are derived via the installed nisarhdf toolkit.


N.Steiner, 2025

"""
from __future__ import annotations

import argparse
import json
import logging
import pathlib
import sys
import time
from typing import Dict, Iterable, Optional, Tuple

import h5py
import numpy as np
import pystac
from osgeo import gdal, gdal_array, ogr, osr
from pystac.extensions.projection import ProjectionExtension
from pystac.utils import str_to_datetime

import builtins  # noqa: E402
import os as _os  # noqa: E402

# nisarhdf modules expect 'os' to be available before their own import cycle.
builtins.os = _os  # type: ignore[attr-defined]

import nisarhdf as nisarhdf_pkg  # noqa: E402
from nisarhdf.nisarGCOVHDF import nisarGCOVHDF  # noqa: E402

# Resolve the actual repo path from the installed package (editable installs work too)
NISARHDF_REPO = str(pathlib.Path(nisarhdf_pkg.__file__).parent)

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS")

PROCESSING_SOFTWARE = "nisarhdf"

# Use a module-level logger that propagates to parent (e.g., Prefect's logger)
logger = logging.getLogger(__name__)

COG_CREATION_OPTIONS = [
    "COMPRESS=DEFLATE",
    "LEVEL=6",
    "BLOCKSIZE=512",
    "OVERVIEWS=IGNORE_EXISTING",
    "SPARSE_OK=YES",
    "NUM_THREADS=ALL_CPUS",
    "RESAMPLING=AVERAGE",
    "PREDICTOR=YES",
    "BIGTIFF=IF_SAFER",
]

COG_TRANSLATE_OPTIONS = gdal.TranslateOptions(format="COG", creationOptions=COG_CREATION_OPTIONS)


def decode_attr(value, *, decode_bytes: bool = False):
    """Normalize HDF5 attribute scalars."""
    if value is None:
        return None
    if isinstance(value, h5py.Dataset):
        value = value[()]
    if isinstance(value, (bytes, np.bytes_)):
        raw = bytes(value)
        return raw.decode("utf-8") if decode_bytes else raw
    if hasattr(value, "item"):
        value = value.item()
    if decode_bytes and isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def ensure_dir(path: pathlib.Path) -> None:
    """Ensure a directory exists or create it."""
    path.mkdir(parents=True, exist_ok=True)


def ensure_band_loaded(reader: nisarGCOVHDF, band: str):
    """Guarantee that a band is resident in memory."""
    data = getattr(reader, band, None)
    if data is None:
        reader.loadData([band], resetFields=False)
        data = getattr(reader, band)
    return data


def derive_units(reader: nisarGCOVHDF, band: str, attr_units: Optional[str]) -> Optional[str]:
    """Mirror nisarGCOVHDF conversions when reporting per-band units."""
    cov_terms = getattr(reader, "covTerms", [])
    if band in cov_terms:
        if getattr(reader, "units", None) == "dB":
            return "dB"
        backscatter = getattr(reader, "backscatterType", None)
        if backscatter in ("sigma0", "gamma0"):
            return backscatter
        return getattr(reader, "units", attr_units)
    return attr_units


def ensure_incidence_angle(reader: nisarGCOVHDF) -> None:
    """Populate incidence angle data and register it for export."""
    if getattr(reader, "incidenceAngle", None) is not None:
        logger.debug("incidenceAngle already present: shape=%s", getattr(reader.incidenceAngle, "shape", None))
        return

    if not hasattr(reader, "xGrid") or not hasattr(reader, "yGrid"):
        logger.debug("Setting up XY grid for incidence angle computation")
        reader.setupXYGrid()

    mask_array = None
    try:
        mask_array = ensure_band_loaded(reader, "mask")
    except Exception:
        mask_array = None

    if mask_array is not None:
        logger.debug("Using mask for noDataLocations (mask==0)")
        reader.noDataLocations = mask_array == 0
    else:
        reader.noDataLocations = None

    t0 = time.monotonic()
    z = np.zeros_like(reader.xGrid, dtype=np.float32)
    incidence = reader.incidenceAngleCube(
        reader.xGrid,
        reader.yGrid,
        z,
        maskNoData=reader.noDataLocations,
        save=False,
    )
    if incidence is None:
        logger.debug("incidenceAngleCube returned None; skipping incidence angle export")
        return

    incidence = np.asarray(incidence, dtype=np.float32)
    if reader.noDataLocations is not None:
        incidence[reader.noDataLocations] = np.nan

    reader.incidenceAngle = incidence
    if "incidenceAngle" not in reader.dataFields:
        reader.dataFields.append("incidenceAngle")
    reader.noDataValuesTiff["incidenceAngle"] = np.nan
    dt = time.monotonic() - t0
    logger.debug(
        "Computed incidenceAngle: shape=%s dtype=%s elapsed=%.3fs",
        incidence.shape,
        incidence.dtype,
        dt,
    )


def apply_metadata(
    ds: gdal.Dataset,
    *args,
    frequency: str,
    subdataset: str,
    epsg_value: Optional[int],
    band_name: str,
    long_name: Optional[str],
    description: Optional[str],
    units: Optional[str],
    no_data: Optional[float],
    backscatter: Optional[str],
     **kwargs
) -> None:
    """Apply dataset and band metadata while the dataset is open."""
    logger.debug(
        "Applying metadata: freq=%s band=%s epsg=%s units=%s nodata=%s backscatter=%s",
        frequency,
        band_name,
        epsg_value,
        units,
        no_data,
        backscatter,
    )
    ds.SetMetadataItem("frequency", frequency)
    ds.SetMetadataItem("source_subdataset", subdataset)
    ds.SetMetadataItem("processing_software", PROCESSING_SOFTWARE)
    ds.SetMetadataItem("processing_software_path", NISARHDF_REPO)
    if backscatter:
        ds.SetMetadataItem("backscatter_type", backscatter)
    if epsg_value is not None:
        ds.SetMetadataItem("proj_epsg", str(epsg_value))
    
    # Apply additional metadata from args and kwargs, to ds level only
    for arg_name, arg_value in args:
        if arg_value is not None:
            ds.SetMetadataItem(arg_name, str(arg_value))
    for arg_name, arg_value in kwargs.items():
        if arg_value is not None:
            ds.SetMetadataItem(arg_name, str(arg_value))

    band = ds.GetRasterBand(1)
    band.SetDescription(long_name or band_name)
    band.SetMetadataItem("band_name", band_name)
    if long_name:
        band.SetMetadataItem("long_name", long_name)
    if description:
        band.SetMetadataItem("description", description)
    if units:
        band.SetMetadataItem("units", units)
    if no_data is not None and not np.isnan(no_data):
        band.SetMetadataItem("nodata_value", str(no_data))

    band.FlushCache()
# what would be the direction of the flight
derived_data = [
   "incidenceAngle", "elevationAngle", "azimuthAngle", "depressionAngle"
]


def export_frequency_assets(
    reader: nisarGCOVHDF,
    freq_name: str,
    h5_path: pathlib.Path,
    output_root: pathlib.Path,
    *,
    overwrite: bool,
    sigma0: bool,
    db: bool,
    skip_incidence_angle: bool = False,
    verbose: bool = False,
    subset_window: Optional[Tuple[int, int, int, int]] = None,
) -> Dict[str, Dict[str, Optional[str]]]:
    """
    Export all bands for a frequency to COG files and return STAC asset metadata.
    
    Data flow:
    - Converted bands (sigma0/dB): already in memory from openHDF, written directly via write_mem_cog
    - Unconverted bands: GDAL reads directly from HDF5 via translate_subdataset_to_cog (no double-read)
    - Derived bands (incidenceAngle): computed on-demand, written via write_mem_cog
    - Fallback: if GDAL fails, load to memory via ensure_band_loaded and write_mem_cog
    """
    

    logger.debug(
        "Exporting frequency assets: freq=%s product=%s bands_group=%s",
        freq_name,
        getattr(reader, "product", None),
        getattr(reader, "bands", None),
    )
    group = reader.h5[reader.product][reader.bands][freq_name]
    freq_dir = output_root / h5_path.stem / freq_name
    ensure_dir(freq_dir)

    if not skip_incidence_angle:
        ensure_incidence_angle(reader)

    epsg_value = int(reader.epsg) if reader.epsg is not None else None
    base_gt = reader.getGeoTransform(tiff=True)
    geo_transform = base_gt
    if subset_window is not None:
        xoff, yoff, xsize, ysize = subset_window
        # Adjust geotransform for window: GT' = GT + [xoff*GT[1] + yoff*GT[2], 0, 0, xoff*GT[4] + yoff*GT[5]] on origins
        geo_transform = (
            float(base_gt[0]) + xoff * float(base_gt[1]) + yoff * float(base_gt[2]),
            float(base_gt[1]),
            float(base_gt[2]),
            float(base_gt[3]) + xoff * float(base_gt[4]) + yoff * float(base_gt[5]),
            float(base_gt[4]),
            float(base_gt[5]),
        )
        logger.debug(
            "Subset window: xoff=%d yoff=%d xsize=%d ysize=%d | adjusted_gt=%s",
            xoff,
            yoff,
            xsize,
            ysize,
            geo_transform,
        )
    logger.debug("Geo: epsg=%s geo_transform=%s", epsg_value, tuple(float(v) for v in geo_transform))
    
    # Prepare derived sources like incidence angle
    derived_sources: Dict[str, np.ndarray] = {}
    derived_attrs: Dict[str, h5py.Dataset] = {}
    if hasattr(reader, "incidenceAngle"):
        derived_sources["incidenceAngle"] = reader.incidenceAngle
        derived_attrs["incidenceAngle"] = reader.h5[reader.product]["metadata"]["radarGrid"]["incidenceAngle"]

    assets: Dict[str, Dict[str, Optional[str]]] = {}
    bands_to_process = list(dict.fromkeys(list(reader.dataFields) + list(derived_sources.keys())))
    cov_terms = set(getattr(reader, "covTerms", []))
    logger.debug(
        "Bands to process (%d): %s | cov_terms=%s",
        len(bands_to_process),
        ",".join(bands_to_process),
        list(cov_terms),
    )

    for band in bands_to_process:
        dataset = group.get(band)
        is_derived = band in derived_sources
        # If conversions were requested, export covariance terms from memory so the
        # converted values (sigma0 and/or dB) are written to the COG instead of raw HDF5.
        force_from_memory = (sigma0 or db) and (band in cov_terms)
        logger.debug(
            "Band=%s derived=%s force_from_memory=%s",
            band,
            is_derived,
            force_from_memory,
        )
        if not is_derived and not force_from_memory:
            if not isinstance(dataset, h5py.Dataset) or dataset.ndim != 2:
                logger.debug("Skipping band=%s (not 2D dataset)", band)
                continue

        dst_path = freq_dir / f"{h5_path.stem}.{band}.tif"
        if overwrite and dst_path.exists():
            dst_path.unlink()

        write_nodata: Optional[float] = None
        if is_derived or force_from_memory:
            # If forcing from memory (converted covariance), pull data from the reader
            # and keep attribute source from the original dataset for metadata.
            if force_from_memory:
                data_array = getattr(reader, band)
                attr_source = dataset  # retain HDF5 attrs like long_name/description
            else:
                data_array = derived_sources[band]
                attr_source = derived_attrs[band]
            # Apply subset window slicing if requested
            if subset_window is not None and data_array is not None:
                xoff, yoff, xsize, ysize = subset_window
                data_array = np.asarray(data_array)[yoff : yoff + ysize, xoff : xoff + xsize]
            write_nodata = reader.findNoDataValue(band, tiff=True)
            dtype_str = str(np.asarray(data_array).dtype)
            subdataset = f"DERIVED:{band}"
            logger.debug(
                "Band=%s write from memory: shape=%s dtype=%s nodata=%s subdataset=%s",
                band,
                np.asarray(data_array).shape,
                dtype_str,
                write_nodata,
                subdataset,
            )
            # Fallback conversions in case upstream didn't apply
            override_units: Optional[str] = None
            override_backscatter: Optional[str] = None
            if band in cov_terms:
                if sigma0 and getattr(reader, "backscatterType", None) != "sigma0":
                    try:
                        rtc = getattr(reader, "rtcGammaToSigmaFactor")
                        if rtc is not None:
                            data_array = np.asarray(data_array, dtype=np.float32) * np.asarray(rtc, dtype=np.float32)
                            override_backscatter = "sigma0"
                    except Exception:
                        pass
                if db and getattr(reader, "units", None) != "dB":
                    with np.errstate(divide="ignore", invalid="ignore"):
                        data_array = 10.0 * np.log10(np.asarray(data_array, dtype=np.float32))
                    override_units = "dB"
            if verbose and (sigma0 or db) and band in cov_terms:
                try:
                    vmin = float(np.nanmin(data_array))
                    vmax = float(np.nanmax(data_array))
                    logger.debug(
                        "writing converted %s from memory: units=%s, backscatter=%s, min/max=%.2f/%.2f",
                        band,
                        getattr(reader, "units", None),
                        getattr(reader, "backscatterType", None),
                        vmin,
                        vmax,
                    )
                except Exception:
                    pass
        else:
            data_array = None
            attr_source = dataset
            no_data_value = reader.findNoDataValue(band, tiff=True)
            write_nodata = no_data_value
            dtype_str = str(dataset.dtype)
            subdataset = (
                f'HDF5:"{h5_path}"://science/LSAR/{reader.product}/grids/{freq_name}/{band}'
            )
            logger.debug(
                "Band=%s stream via GDAL: dtype=%s nodata=%s subdataset=%s",
                band,
                dtype_str,
                write_nodata,
                subdataset,
            )

        attr_long = decode_attr(attr_source.attrs.get("long_name"), decode_bytes=True) if attr_source is not None else None
        attr_desc = decode_attr(attr_source.attrs.get("description"), decode_bytes=True) if attr_source is not None else None
        attr_units = decode_attr(attr_source.attrs.get("units"), decode_bytes=True) if attr_source is not None else None
        units = derive_units(reader, band, attr_units)
        backscatter = getattr(reader, "backscatterType", None) if band in getattr(reader, "covTerms", []) else None
        # Apply any overrides from fallback conversions
        if 'override_units' in locals() and override_units is not None:
            units = override_units
        if 'override_backscatter' in locals() and override_backscatter is not None:
            backscatter = override_backscatter

        stac_nodata = None
        if write_nodata is not None and not np.isnan(float(write_nodata)):
            stac_nodata = float(write_nodata)

        should_write = overwrite or not dst_path.exists()

        logger.debug(
            "Output target: %s (should_write=%s) units=%s backscatter=%s",
            dst_path,
            should_write,
            units,
            backscatter,
        )
        if should_write:
            if data_array is not None:
                logger.debug("Writing band=%s from memory to %s", band, dst_path)
                write_mem_cog(
                    dst_path,
                    data_array,
                    geo_transform,
                    epsg_value,
                    write_nodata,
                    frequency=freq_name,
                    subdataset=subdataset,
                    band_name=band,
                    long_name=attr_long,
                    description=attr_desc,
                    units=units,
                    backscatter=backscatter,
                )
            else:
                try:
                    logger.debug("Translating subdataset to COG: %s -> %s", subdataset, dst_path)
                    translate_subdataset_to_cog(
                        subdataset,
                        dst_path,
                        geo_transform,
                        epsg_value,
                        write_nodata,
                        frequency=freq_name,
                        band_name=band,
                        long_name=attr_long,
                        description=attr_desc,
                        units=units,
                        backscatter=backscatter,
                        srcwin=subset_window,
                    )
                except RuntimeError as err:
                    logger.warning(
                        "%s/%s: %s — falling back to in-memory write",
                        freq_name,
                        band,
                        err,
                    )
                    data_array = ensure_band_loaded(reader, band)
                    dtype_str = str(data_array.dtype)
                    logger.debug("Fallback load band=%s dtype=%s shape=%s", band, dtype_str, np.asarray(data_array).shape)
                    write_mem_cog(
                        dst_path,
                        data_array,
                        geo_transform,
                        epsg_value,
                        write_nodata,
                        frequency=freq_name,
                        subdataset=subdataset,
                        band_name=band,
                        long_name=attr_long,
                        description=attr_desc,
                        units=units,
                        backscatter=backscatter,
                    )
        else:
            logger.info("Skipping existing %s", dst_path)

        assets[f"{freq_name}_{band}"] = {
            "href": str(dst_path.resolve()),
            "title": attr_long or band,
            "description": attr_desc,
            "units": units,
            "dtype": dtype_str,
            "no_data": stac_nodata,
            "frequency": freq_name,
            "proj:epsg": epsg_value,
            "processing:software": PROCESSING_SOFTWARE,
            "processing:backscatter": backscatter,
        }

    return assets


def write_mem_cog(
    dst_path: pathlib.Path,
    data: np.ndarray,
    geo_transform,
    epsg_value: Optional[int],
    no_data_value: Optional[float],
    *,
    frequency: str,
    subdataset: str,
    band_name: str,
    long_name: Optional[str],
    description: Optional[str],
    units: Optional[str],
    backscatter: Optional[str],
) -> None:
    """Write a numpy array to COG using GDAL in-memory dataset."""
    logger.debug(
        "write_mem_cog: dst=%s shape=%s dtype=%s epsg=%s nodata=%s",
        dst_path,
        np.asarray(data).shape,
        np.asarray(data).dtype,
        epsg_value,
        no_data_value,
    )
    mem_driver = gdal.GetDriverByName("MEM")
    array = np.asarray(data)
    rows, cols = array.shape
    gdal_type = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)
    mem_ds = mem_driver.Create("", cols, rows, 1, gdal_type)
    mem_ds.SetGeoTransform(tuple(float(v) for v in geo_transform))
    if epsg_value is not None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg_value)
        mem_ds.SetProjection(srs.ExportToWkt())
    mem_band = mem_ds.GetRasterBand(1)
    mem_band.WriteArray(array)
    if no_data_value is not None:
        mem_band.SetNoDataValue(float(no_data_value))
    mem_band.FlushCache()

    out_ds = gdal.Translate(str(dst_path), mem_ds, options=COG_TRANSLATE_OPTIONS)
    if out_ds is None:
        raise RuntimeError(f"GDAL translation failed for {dst_path}")
    out_ds.SetGeoTransform(tuple(float(v) for v in geo_transform))
    if epsg_value is not None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg_value)
        out_ds.SetProjection(srs.ExportToWkt())
    if no_data_value is not None:
        out_ds.GetRasterBand(1).SetNoDataValue(float(no_data_value))
    apply_metadata(
        out_ds,
        frequency=frequency,
        subdataset=subdataset,
        epsg_value=epsg_value,
        band_name=band_name,
        long_name=long_name,
        description=description,
        units=units,
        no_data=no_data_value,
        backscatter=backscatter,
    )
    out_ds.FlushCache()
    logger.debug("write_mem_cog complete: %s", dst_path)
    out_ds = None
    mem_band = None
    mem_ds = None


def translate_subdataset_to_cog(
    subdataset: str,
    dst_path: pathlib.Path,
    geo_transform,
    epsg_value: Optional[int],
    no_data_value: Optional[float],
    *,
    frequency: str,
    band_name: str,
    long_name: Optional[str],
    description: Optional[str],
    units: Optional[str],
    backscatter: Optional[str],
    srcwin: Optional[Tuple[int, int, int, int]] = None,
) -> None:
    """Use GDAL to translate an HDF5 subdataset directly to COG if supported."""
    logger.debug("translate_subdataset_to_cog: %s -> %s", subdataset, dst_path)
    src = gdal.Open(subdataset, gdal.GA_ReadOnly)
    if src is None:
        raise RuntimeError(f"GDAL could not open {subdataset}")

    if srcwin is not None:
        xoff, yoff, xsize, ysize = srcwin
        logger.debug("GDAL translate with srcWin=%s", (xoff, yoff, xsize, ysize))
        translate_opts = gdal.TranslateOptions(
            format="COG", creationOptions=COG_CREATION_OPTIONS, srcWin=[xoff, yoff, xsize, ysize]
        )
    else:
        translate_opts = COG_TRANSLATE_OPTIONS

    out_ds = gdal.Translate(str(dst_path), src, options=translate_opts)
    if out_ds is None:
        raise RuntimeError(f"GDAL translation failed for {dst_path}")

    # Ensure georeferencing is present even if the driver omitted it.
    out_ds.SetGeoTransform(tuple(float(v) for v in geo_transform))
    if epsg_value is not None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg_value)
        out_ds.SetProjection(srs.ExportToWkt())
    if no_data_value is not None:
        out_ds.GetRasterBand(1).SetNoDataValue(float(no_data_value))

    apply_metadata(
        out_ds,
        frequency=frequency,
        subdataset=subdataset,
        epsg_value=epsg_value,
        band_name=band_name,
        long_name=long_name,
        description=description,
        units=units,
        no_data=no_data_value,
        backscatter=backscatter,
    )

    out_ds.FlushCache()
    logger.debug("translate_subdataset_to_cog complete: %s", dst_path)
    out_ds = None
    src = None


def geometry_from_wkt(
    wkt: Optional[str], fallback_bbox: Tuple[float, float, float, float]
) -> Tuple[Dict, list]:
    if wkt:
        geom = ogr.CreateGeometryFromWkt(wkt)
        if geom:
            geom.FlattenTo2D()
            envelope = geom.GetEnvelope()  # minx, maxx, miny, maxy
            bbox = [envelope[0], envelope[2], envelope[1], envelope[3]]
            geom_dict = json.loads(geom.ExportToJson())
            return geom_dict, bbox
    minx, miny, maxx, maxy = fallback_bbox
    geometry = {
        "type": "Polygon",
        "coordinates": [
            [
                [minx, miny],
                [minx, maxy],
                [maxx, maxy],
                [maxx, miny],
                [minx, miny],
            ]
        ],
    }
    return geometry, [minx, miny, maxx, maxy]


def fallback_bbox_from_reader(reader: nisarGCOVHDF) -> Tuple[float, float, float, float]:
    xs = np.asarray(reader.xCoords)
    ys = np.asarray(reader.yCoords)
    return float(xs.min()), float(ys.min()), float(xs.max()), float(ys.max())


def extract_item_metadata(reader: nisarGCOVHDF) -> Dict[str, Optional[str]]:
    ident = reader.h5["identification"]
    item_id = decode_attr(ident["granuleId"][()], decode_bytes=True)
    start = decode_attr(ident["zeroDopplerStartTime"][()], decode_bytes=True)
    end = decode_attr(ident["zeroDopplerEndTime"][()], decode_bytes=True)
    mission = decode_attr(ident.get("missionId"), decode_bytes=True)
    product = decode_attr(ident.get("productType"), decode_bytes=True)
    wkt = decode_attr(ident.get("boundingPolygon"), decode_bytes=True)
    fallback_bbox = fallback_bbox_from_reader(reader)
    geometry, bbox = geometry_from_wkt(wkt, fallback_bbox)

    start_dt = str_to_datetime(start) if start else None
    end_dt = str_to_datetime(end) if end else None
    epsg_value = int(reader.epsg) if reader.epsg is not None else None

    return {
        "id": item_id,
        "geometry": geometry,
        "bbox": bbox,
        "start": start_dt,
        "end": end_dt,
        "epsg": epsg_value,
        "properties": {
            "mission": mission,
            "product": product,
            "processing:software": PROCESSING_SOFTWARE,
            "processing:repository": NISARHDF_REPO,
        },
    }


def open_frequency_readers(
    h5_path: pathlib.Path,
    *,
    sigma0: bool = False,
    db: bool = False,
    downsample_factor: int = 1,
) -> Dict[str, nisarGCOVHDF]:
    """
    Open GCOV readers for all frequencies in the HDF5 file.
    
    Data loading strategy:
    - If sigma0 or db conversions are requested, data is loaded eagerly (noLoadData=False)
      because nisarGCOVHDF applies conversions during openHDF when those flags are set.
    - Otherwise (noLoadData=True), data stays on disk and is loaded lazily when first accessed.
    
    All readers use consistent downsample_factor to avoid dimension mismatches.
    """
    logger.debug(
        "open_frequency_readers: path=%s sigma0=%s db=%s downsample=%s",
        h5_path,
        sigma0,
        db,
        downsample_factor,
    )
    readers: Dict[str, nisarGCOVHDF] = {}
    no_load_data = not (sigma0 or db)
    base_reader = nisarGCOVHDF()
    downsample_factor_inp = {'downsampleFactorRow': downsample_factor, 'downsampleFactorColumn': downsample_factor}
    base_reader.openHDF(
        str(h5_path),
        noLoadData=no_load_data, sigma0=sigma0, dB=db, downsampleFactor=downsample_factor_inp
    )
    readers[base_reader.frequency] = base_reader
    grid_group = base_reader.h5[base_reader.product]["grids"]
    logger.debug("Base frequency=%s, available=%s", base_reader.frequency, list(grid_group.keys()))
    for freq in list(grid_group.keys()):
        if freq in readers:
            continue
        freq_reader = nisarGCOVHDF(frequency=freq)
        freq_reader.openHDF(
            str(h5_path),
            noLoadData=no_load_data,
            sigma0=sigma0,
            dB=db,
            downsampleFactor=downsample_factor_inp
        )
        readers[freq] = freq_reader
    logger.debug("Opened %d readers: %s", len(readers), list(readers.keys()))
    return readers


def close_readers(readers: Dict[str, nisarGCOVHDF]) -> None:
    for reader in set(readers.values()):
        try:
            reader.close()
        except Exception:
            pass


def build_item(
    metadata: Optional[Dict[str, Optional[str]]],
    assets: Dict[str, Dict[str, Optional[str]]],
    stac_dir: pathlib.Path,
) -> None:
    if not metadata or not assets:
        logger.debug("build_item: no metadata or assets; skipping STAC creation")
        return

    item_datetime = metadata.get("start") or metadata.get("end")
    item = pystac.Item(
        id=metadata["id"],
        geometry=metadata["geometry"],
        bbox=metadata["bbox"],
        datetime=item_datetime,
        properties={k: v for k, v in metadata["properties"].items() if v is not None},
    )
    if metadata["start"]:
        item.common_metadata.start_datetime = metadata["start"]
    if metadata["end"]:
        item.common_metadata.end_datetime = metadata["end"]
    if metadata.get("epsg") is not None:
        ProjectionExtension.ext(item, add_if_missing=True).epsg = metadata["epsg"]

    for key, info in assets.items():
        asset = pystac.Asset(
            href=info["href"],
            media_type=pystac.MediaType.COG,
            roles=["data"],
            title=info.get("title"),
            description=info.get("description"),
        )
        extra_fields = {
            "frequency": info.get("frequency"),
            "raster:units": info.get("units"),
            "raster:datatype": info.get("dtype"),
            "proj:epsg": info.get("proj:epsg"),
            "processing:software": info.get("processing:software"),
            "processing:backscatter": info.get("processing:backscatter"),
        }
        if info.get("no_data") is not None:
            extra_fields["raster:nodata"] = info["no_data"]
        asset.extra_fields.update({k: v for k, v in extra_fields.items() if v is not None})
        item.add_asset(key, asset)

    ensure_dir(stac_dir)
    item_path = stac_dir / f"{metadata['id']}.json"
    item.set_self_href(str(item_path))
    item.save_object(include_self_link=False, dest_href=str(item_path))
    logger.debug("STAC item written: %s (assets=%d)", item_path, len(assets))


def parse_args(argv: Optional[Iterable[str]] = None):
    parser = argparse.ArgumentParser(
        description="Export NISAR L2 GCOV bands to Cloud-Optimized GeoTIFF and generate PySTAC items."
    )
    parser.add_argument("inputs", nargs="+", type=pathlib.Path, help="Paths to GCOV .h5 granules.")
    
    
    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("gcov_cogs"),
        help="Destination directory for COGs (default: ./gcov_cogs).",
    )
    parser.add_argument(
        "--stac",
        type=pathlib.Path,
        default=pathlib.Path("stac_items"),
        help="Output directory for PySTAC item JSON (default: ./stac_items).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs.",
    )
    parser.add_argument(
        "--sigma0",
        action="store_true",
        help="Convert gamma0 covariance terms to sigma0 (matches nisarh5toimage GCOV option).",
    )
    parser.add_argument(
        "--db",
        action="store_true",
        help="Convert backscatter covariance terms to decibels after other conversions.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose debug logging (prints band-level conversion info).",
    )
    parser.add_argument(
        "--no-incidence-angle",
        action="store_true",
        help="Skip incidence angle computation and export (faster processing).",
    )
    parser.add_argument(
        "--downsample-factor",
        type=int,
        default=1,
        help="Downsample factor for output COGs (default: 1, no downsampling).",
    )
    parser.add_argument(
        "--srcwin",
        nargs=4,
        type=int,
        metavar=("XOFF", "YOFF", "XSIZE", "YSIZE"),
        help="Subset read window in pixels: x offset, y offset, width, height.",
    )

    args = parser.parse_args(argv)
    if getattr(args, "srcwin", None) is not None:
        # Convert list to tuple for downstream use
        args.srcwin = tuple(int(v) for v in args.srcwin)  # type: ignore[attr-defined]
    # Note: logging may not yet be configured; defer rich debug to main after setup
    return args


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(argv)
    ensure_dir(args.output)
    ensure_dir(args.stac)
    logger.debug(
        "Args: inputs=%s output=%s stac=%s overwrite=%s sigma0=%s db=%s verbose=%s no_incidence=%s downsample=%d srcwin=%s",
        [p.name for p in args.inputs],
        args.output,
        args.stac,
        args.overwrite,
        args.sigma0,
        args.db,
        args.verbose,
        args.no_incidence_angle,
        args.downsample_factor,
        getattr(args, "srcwin", None),
    )

    for h5_path in args.inputs:
        if not h5_path.exists():
            logger.error("Input not found: %s", h5_path)
            continue

        logger.info("Processing %s", h5_path.name)
        assets: Dict[str, Dict[str, Optional[str]]] = {}
        readers: Dict[str, nisarGCOVHDF] = {}
        item_metadata: Optional[Dict[str, Optional[str]]] = None

        try:
            readers = open_frequency_readers(
                h5_path, 
                sigma0=args.sigma0, 
                db=args.db, 
                downsample_factor=args.downsample_factor
                )
            base_reader = readers.get("frequencyA") or next(iter(readers.values()))
            item_metadata = extract_item_metadata(base_reader)
            logger.debug(
                "Item metadata extracted: id=%s epsg=%s time=[%s,%s]",
                item_metadata.get("id") if item_metadata else None,
                item_metadata.get("epsg") if item_metadata else None,
                item_metadata.get("start") if item_metadata else None,
                item_metadata.get("end") if item_metadata else None,
            )

            for freq_name in readers.keys():
                reader = readers[freq_name]
                try:
                    freq_assets = export_frequency_assets(
                        reader,
                        freq_name,
                        h5_path,
                        args.output,
                        overwrite=args.overwrite,
                        sigma0=args.sigma0,
                        db=args.db,
                        skip_incidence_angle=args.no_incidence_angle,
                        verbose=args.verbose,
                        subset_window=getattr(args, "srcwin", None),
                    )
                except Exception as exc:
                    logger.warning("%s: %s", freq_name, exc)
                    continue
                assets.update(freq_assets)
            logger.debug("Finished %s: wrote %d assets", h5_path.name, len(assets))
        except Exception as exc:
            logger.error("%s: %s", h5_path.name, exc)
            assets = {}
            item_metadata = None
        finally:
            close_readers(readers)

        build_item(item_metadata, assets, args.stac)

    return 0


def setup_logging(verbose: bool = False) -> None:
    """
    Configure logging for standalone CLI usage with rich pretty formatting.
    When invoked from Prefect (or any framework), skip this and let the parent configure.
    """
    if logging.getLogger().hasHandlers():
        # Already configured (e.g., by Prefect), respect parent setup
        return
    
    level = logging.DEBUG if verbose else logging.INFO
    
    try:
        from rich.logging import RichHandler
        
        logging.basicConfig(
            level=level,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(
                rich_tracebacks=True,
                tracebacks_show_locals=verbose,
                markup=True,
                show_time=True,
                show_level=True,
                show_path=verbose,
            )]
        )
    except ImportError:
        # Fallback to basic logging if rich is not available
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )


if __name__ == "__main__":
    # Parse args early to set up logging before main() runs
    import sys
    args = parse_args(sys.argv[1:])
    setup_logging(verbose=args.verbose)
    raise SystemExit(main(sys.argv[1:]))
