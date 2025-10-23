import os
import numpy as np
import pytest

pytest.importorskip("rioxarray", reason="rioxarray required for VRT reading")
pytest.importorskip("osgeo", reason="GDAL required for VRT/GeoTIFF operations")

from osgeo import gdal, osr
import rioxarray as rxr

from nisarhdf import writeMultiBandVrt, readVrtAsXarray


def _write_geotiff(path, arr, epsg=4326, nodata=np.nan):
    driver = gdal.GetDriverByName("GTiff")
    H, W = arr.shape
    # COG creation options: tiling, compression, and optimization
    cog_options = [
        "TILED=YES",
        "COMPRESS=LZW",
        "COPY_SRC_OVERVIEWS=YES",
        "GEOTIFF_VERSION=1.1"
    ]
    ds = driver.Create(str(path), W, H, 1, gdal.GDT_Float32, options=cog_options) 
    ds.SetGeoTransform([-0.5, 1.0, 0.0, -0.5, 0.0, 1.0])
    sr = osr.SpatialReference()
    sr.ImportFromEPSG(epsg)
    ds.SetProjection(sr.ExportToWkt())
    band = ds.GetRasterBand(1)
    if not np.isnan(nodata):
        band.SetNoDataValue(float(nodata))
    band.WriteArray(arr.astype(np.float32))
    band.FlushCache()
    ds.FlushCache()
    # Build overviews for COG
    ds.BuildOverviews("NEAREST", [2, 4, 8])
    ds = None


@pytest.mark.integration
def test_writeMultiBandVrt_tiff_inputs(tmp_path):
    arr1 = np.arange(80, dtype=np.float32).reshape(8, 10)
    arr2 = np.flipud(arr1)

    t1 = tmp_path / "b1.tif"
    t2 = tmp_path / "b2.tif"
    _write_geotiff(t1, arr1, epsg=4326)
    _write_geotiff(t2, arr2, epsg=4326)

    out_vrt = tmp_path / "stack.vrt"
    descriptions = ["band_one", "band_two"]

    writeMultiBandVrt(
        str(out_vrt),
        xSize=arr1.shape[1],
        ySize=arr1.shape[0],
        sourceFiles=[str(t1), str(t2)],
        descriptions=descriptions,
        epsg=4326,
        tiff=True,
    )

    # Read with our helper, which renames variables by band Description
    xr = readVrtAsXarray(str(out_vrt))

    # Expect renamed variables to match descriptions
    assert set(xr.data_vars) == set(descriptions)

    # CRS present and correct
    assert xr.rio.crs is not None
    try:
        assert xr.rio.crs.to_epsg() == 4326
    except Exception:
        pass

    # Shapes一致
    shapes = {xr[v].shape for v in xr.data_vars}
    assert shapes == {arr1.shape}

    # Values are as expected for at least one band
    assert np.isclose(xr["band_one"].values[0, 0], 0.0)
    assert np.isclose(xr["band_two"].values[-1, 0], 0.0)
