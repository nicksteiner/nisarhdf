#!/usr/bin/env python3
"""
Convert NISAR L2 GCOV granules into Cloud-Optimized GeoTIFFs and PySTAC items.
Georeferencing and metadata are derived via the nisarhdf toolkit located at
/home/ubuntu/nisarhdf to ensure consistency with mission utilities.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Dict, Iterable, Optional, Tuple

import h5py
import numpy as np
import pystac
from osgeo import gdal, gdal_array, ogr, osr
from pystac.extensions.projection import ProjectionExtension
from pystac.utils import str_to_datetime

NISARHDF_REPO = "/home/ubuntu/nisarhdf"
if NISARHDF_REPO not in sys.path:
    sys.path.insert(0, NISARHDF_REPO)

import builtins  # noqa: E402
import os as _os  # noqa: E402

# nisarhdf modules expect 'os' to be available before their own import cycle.
builtins.os = _os  # type: ignore[attr-defined]

from nisarhdf.nisarGCOVHDF import nisarGCOVHDF  # noqa: E402

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_NUM_THREADS", "ALL_CPUS")

PROCESSING_SOFTWARE = "nisarhdf"

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
        return

    if not hasattr(reader, "xGrid") or not hasattr(reader, "yGrid"):
        reader.setupXYGrid()

    mask_array = None
    try:
        mask_array = ensure_band_loaded(reader, "mask")
    except Exception:
        mask_array = None

    if mask_array is not None:
        reader.noDataLocations = mask_array == 0
    else:
        reader.noDataLocations = None

    z = np.zeros_like(reader.xGrid, dtype=np.float32)
    incidence = reader.incidenceAngleCube(
        reader.xGrid,
        reader.yGrid,
        z,
        maskNoData=reader.noDataLocations,
        save=False,
    )
    if incidence is None:
        return

    incidence = np.asarray(incidence, dtype=np.float32)
    if reader.noDataLocations is not None:
        incidence[reader.noDataLocations] = np.nan

    reader.incidenceAngle = incidence
    if "incidenceAngle" not in reader.dataFields:
        reader.dataFields.append("incidenceAngle")
    reader.noDataValuesTiff["incidenceAngle"] = np.nan


def apply_metadata(
    ds: gdal.Dataset,
    *,
    frequency: str,
    subdataset: str,
    epsg_value: Optional[int],
    band_name: str,
    long_name: Optional[str],
    description: Optional[str],
    units: Optional[str],
    no_data: Optional[float],
    backscatter: Optional[str],
) -> None:
    """Apply dataset and band metadata while the dataset is open."""
    ds.SetMetadataItem("frequency", frequency)
    ds.SetMetadataItem("source_subdataset", subdataset)
    ds.SetMetadataItem("processing_software", PROCESSING_SOFTWARE)
    ds.SetMetadataItem("processing_software_path", NISARHDF_REPO)
    if backscatter:
        ds.SetMetadataItem("backscatter_type", backscatter)
    if epsg_value is not None:
        ds.SetMetadataItem("proj_epsg", str(epsg_value))

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


def export_frequency_assets(
    reader: nisarGCOVHDF,
    freq_name: str,
    h5_path: pathlib.Path,
    output_root: pathlib.Path,
    *,
    overwrite: bool,
    sigma0: bool,
    db: bool,
) -> Dict[str, Dict[str, Optional[str]]]:
    group = reader.h5[reader.product][reader.bands][freq_name]
    freq_dir = output_root / h5_path.stem / freq_name
    ensure_dir(freq_dir)

    ensure_incidence_angle(reader)

    epsg_value = int(reader.epsg) if reader.epsg is not None else None
    geo_transform = reader.getGeoTransform(tiff=True)
    derived_sources: Dict[str, np.ndarray] = {}
    derived_attrs: Dict[str, h5py.Dataset] = {}
    if hasattr(reader, "incidenceAngle"):
        derived_sources["incidenceAngle"] = reader.incidenceAngle
        derived_attrs["incidenceAngle"] = reader.h5[reader.product]["metadata"]["radarGrid"]["incidenceAngle"]

    assets: Dict[str, Dict[str, Optional[str]]] = {}
    bands_to_process = list(dict.fromkeys(list(reader.dataFields) + list(derived_sources.keys())))

    for band in bands_to_process:
        dataset = group.get(band)
        is_derived = band in derived_sources
        if not is_derived:
            if not isinstance(dataset, h5py.Dataset) or dataset.ndim != 2:
                continue

        dst_path = freq_dir / f"{h5_path.stem}.{band}.tif"
        if overwrite and dst_path.exists():
            dst_path.unlink()

        write_nodata: Optional[float] = None
        if is_derived:
            data_array = derived_sources[band]
            attr_source = derived_attrs[band]
            write_nodata = np.nan
            dtype_str = str(data_array.dtype)
            subdataset = f"DERIVED:{band}"
        else:
            data_array = None
            attr_source = dataset
            no_data_value = reader.findNoDataValue(band, tiff=True)
            write_nodata = no_data_value
            dtype_str = str(dataset.dtype)
            subdataset = (
                f'HDF5:"{h5_path}"://science/LSAR/{reader.product}/grids/{freq_name}/{band}'
            )

        attr_long = decode_attr(attr_source.attrs.get("long_name"), decode_bytes=True) if attr_source is not None else None
        attr_desc = decode_attr(attr_source.attrs.get("description"), decode_bytes=True) if attr_source is not None else None
        attr_units = decode_attr(attr_source.attrs.get("units"), decode_bytes=True) if attr_source is not None else None
        units = derive_units(reader, band, attr_units)
        backscatter = getattr(reader, "backscatterType", None) if band in getattr(reader, "covTerms", []) else None

        stac_nodata = None
        if write_nodata is not None and not np.isnan(float(write_nodata)):
            stac_nodata = float(write_nodata)

        should_write = overwrite or not dst_path.exists()

        if should_write:
            if data_array is not None:
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
                    )
                except RuntimeError as err:
                    print(
                        f"  [warn] {freq_name}/{band}: {err} — falling back to in-memory write",
                        file=sys.stderr,
                    )
                    data_array = ensure_band_loaded(reader, band)
                    dtype_str = str(data_array.dtype)
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
            print(f"Skipping existing {dst_path}", file=sys.stderr)

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
) -> None:
    """Use GDAL to translate an HDF5 subdataset directly to COG if supported."""
    src = gdal.Open(subdataset, gdal.GA_ReadOnly)
    if src is None:
        raise RuntimeError(f"GDAL could not open {subdataset}")

    out_ds = gdal.Translate(str(dst_path), src, options=COG_TRANSLATE_OPTIONS)
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
) -> Dict[str, nisarGCOVHDF]:
    readers: Dict[str, nisarGCOVHDF] = {}
    no_load_data = not (sigma0 or db)
    base_reader = nisarGCOVHDF()
    base_reader.openHDF(str(h5_path), noLoadData=no_load_data, sigma0=sigma0, dB=db)
    readers[base_reader.frequency] = base_reader
    grid_group = base_reader.h5[base_reader.product]["grids"]
    for freq in list(grid_group.keys()):
        if freq in readers:
            continue
        freq_reader = nisarGCOVHDF(frequency=freq)
        freq_reader.openHDF(str(h5_path), noLoadData=no_load_data, sigma0=sigma0, dB=db)
        readers[freq] = freq_reader
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
    return parser.parse_args(argv)


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = parse_args(argv)
    ensure_dir(args.output)
    ensure_dir(args.stac)

    for h5_path in args.inputs:
        if not h5_path.exists():
            print(f"Input not found: {h5_path}", file=sys.stderr)
            continue

        print(f"Processing {h5_path.name}…")
        assets: Dict[str, Dict[str, Optional[str]]] = {}
        readers: Dict[str, nisarGCOVHDF] = {}
        item_metadata: Optional[Dict[str, Optional[str]]] = None

        try:
            readers = open_frequency_readers(h5_path, sigma0=args.sigma0, db=args.db)
            base_reader = readers.get("frequencyA") or next(iter(readers.values()))
            item_metadata = extract_item_metadata(base_reader)

            for freq_name in sorted(readers.keys()):
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
                    )
                except Exception as exc:
                    print(f"  [warn] {freq_name}: {exc}", file=sys.stderr)
                    continue
                assets.update(freq_assets)
        except Exception as exc:
            print(f"  [error] {h5_path.name}: {exc}", file=sys.stderr)
            assets = {}
            item_metadata = None
        finally:
            close_readers(readers)

        build_item(item_metadata, assets, args.stac)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
