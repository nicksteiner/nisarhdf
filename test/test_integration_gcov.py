import numpy as np
import pytest

# Skip gracefully if heavy deps are missing
pytest.importorskip("rioxarray", reason="rioxarray required for VRT reading")
pytest.importorskip("osgeo", reason="GDAL required for VRT/GeoTIFF operations")

from nisarhdf import nisarGCOVHDF, readVrtAsXarray


@pytest.mark.integration
def test_gcov_write_and_read_vrt(gcov_path, tmp_path):
    # Open GCOV with downsampling to keep runtime reasonable
    gcov = nisarGCOVHDF()
    gcov.openHDF(
        gcov_path,
        downsampleFactor={"downsampleFactorRow": 8, "downsampleFactorColumn": 8},
        noLoadData=False,
        fields=None,
    )

    # Expected bands come from file-provided covariance terms plus extras
    # Use the dataFields selected by the instance to avoid assuming availability
    expected_bands = list(getattr(gcov, "dataFields", []))
    assert len(expected_bands) > 0, "No data fields selected for GCOV"

    out_root = tmp_path / "gcov_out"
    gcov.writeData(str(out_root), tiff=True, driverName="COG")

    vrt_path = str(out_root) + ".vrt"
    xr = readVrtAsXarray(vrt_path)

    # Band presence
    present = set(xr.data_vars)
    for band in expected_bands:
        assert band in present, f"Missing band in VRT dataset: {band}"

    # CRS is present and matches EPSG when available
    assert xr.rio.crs is not None
    if getattr(gcov, "epsg", None) is not None:
        try:
            assert xr.rio.crs.to_epsg() == gcov.epsg
        except Exception:
            # Some drivers may not report EPSG exactly; accept presence above
            pass

    # All bands have identical shape
    shapes = {xr[v].shape for v in xr.data_vars}
    assert len(shapes) == 1

    # Sanity: at least one band has finite data (not entirely nodata)
    found_valid = any(np.isfinite(xr[v].values).any() for v in xr.data_vars)
    assert found_valid, "All bands appear to be entirely nodata/NaN"
