# nisarHDF

This repository contains python code to read the NISAR HDF formats and parse out the necessary parameters and reformatting necessary for the Greenland Ice Mapping Project (*GrIMP*) workflow. It also provide a command line utility (`nisarh5toimage`) for reformatting and extracting NISAR HDF5 data. Most of the functionality though is independent of *GrIMP* and may prove useful for anyone working with NISAR HDF data.

## Python Package
In particular, it simplifies access to NISAR data by allowing:
  - Access data (e.g., `myRUNW.unwrappedPhase`) as numpy variables and other commonly used parameters (e.g., `myGUNW.NumberRangeLooks`) rather than hdf structures,
  - Mask data (e.g., remove  data with connected component values of zero),
  - Apply ionospheric corrections to phase data,
  - Interpolate metadata cubes (e.g., incidenceAngle and baseline), including mapping these variables to the same grid as the data products.
  - Save data and interpolated metadata cubes to binary or geotiff files, and
  - Write a .vrt file to link multiple data layers so they can be read back in as `xarray`.

At present classes have been developed for **RSLC**, **ROFF**, **RUNW**, **RIFG**, **GOFF**, **GUNW**, and **GCOV** products. A tutorial notebook for each of these data types is included in the [Notebook](https://github.com/fastice/nisarhdf/blob/main/Notebooks) folder of this repository.

## Command line utility nisarh5toimage

 While HDF5 is a powerful format for archiving, other formats such as tiff and simple binary can be easier for many users to work with. To facilitate working with NISAR data for these users, the `nisarh5toimage` command line utility can:

- Provide a quick summary of what's in the file (bands etc) along with essential meta data (e.g., *time, orbit, frame, track, center lat/lon*, etc.)
- Pull the individual bands from NISAR hdfs and write them to tiff or binary files.
- Produce an accompanying `.vrt` field that makes it easy to import multiband products into Python (e.g., `rioxarray`) or GIS packages.
- Downsample products (useful for large hi-res GCOV products).
- Produce quickLook PNG files with `.vrt` for geospatial information.


The [nisarh5toimageTutorial](https://github.com/fastice/nisarhdf/blob/main/Notebooks/nisarh5toimageTutorial.ipynb) notebook in the [Notebook](https://github.com/fastice/nisarhdf/blob/main/Notebooks) directory of this respository includes several examples of how to use this program.

## GCOV to COG Conversion Script

The `scripts/nisar_gcov_cog.py` script converts NISAR L2 GCOV HDF5 files to Cloud-Optimized GeoTIFFs (COGs) and generates PySTAC metadata items.

### Example Usage

**Basic conversion (default gamma0 power units):**
```bash
python ./scripts/nisar_gcov_cog.py /path/to/GCOV.h5
```

**Convert to sigma0 backscatter:**
```bash
python ./scripts/nisar_gcov_cog.py --sigma0 /path/to/GCOV.h5
```

**Convert to dB scale:**
```bash
python ./scripts/nisar_gcov_cog.py --db /path/to/GCOV.h5
```

**Combine sigma0 and dB (sigma0 in dB):**
```bash
python ./scripts/nisar_gcov_cog.py --sigma0 --db /path/to/GCOV.h5
```

**Skip expensive incidence angle computation:**
```bash
python ./scripts/nisar_gcov_cog.py --no-incidence-angle /path/to/GCOV.h5
```

**Specify output directory and overwrite existing files:**
```bash
python ./scripts/nisar_gcov_cog.py -o /output/cogs --stac /output/stac --overwrite /path/to/GCOV.h5
```

**Process multiple files with verbose debug output:**
```bash
python ./scripts/nisar_gcov_cog.py --db --no-incidence-angle --verbose -o /tmp/gcov_cogs --overwrite \
  /path/to/GCOV1.h5 \
  /path/to/GCOV2.h5 \
  /path/to/GCOV3.h5
```

**Full workflow example (recommended):**
```bash
# Convert to sigma0 dB, skip incidence angle for speed, enable verbose logging
python ./scripts/nisar_gcov_cog.py \
  --sigma0 \
  --db \
  --no-incidence-angle \
  --verbose \
  -o /data/gcov_cogs \
  --stac /data/stac_items \
  --overwrite \
  /path/to/NISAR_L2_PR_GCOV_*.h5
```

### Output Structure

For an input file `NISAR_L2_PR_GCOV_XXX.h5`, the script produces:

**COG Files** (in `--output` directory):
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.HHHH.tif`
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.VVVV.tif`
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.HVHV.tif`
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.VHVH.tif`
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.HHVV.tif`
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.mask.tif`
- `gcov_cogs/NISAR_L2_PR_GCOV_XXX/frequencyA/NISAR_L2_PR_GCOV_XXX.incidenceAngle.tif` (if not skipped)

**PySTAC Item** (in `--stac` directory):
- `stac_items/NISAR_L2_PR_GCOV_XXX.json`

### Conversion Options

- **Default**: gamma0 backscatter in power (linear) units
- **`--sigma0`**: Convert gamma0 → sigma0 using RTC correction factor
- **`--db`**: Convert power → dB scale (10×log₁₀)
- **`--sigma0 --db`**: Convert gamma0 → sigma0 → dB
- **`--no-incidence-angle`**: Skip incidence angle generation (faster, smaller output)
- **`--verbose`**: Print min/max values and units for each band at write time (debugging)

### Notes

- All outputs are COG-compliant with tiling, compression (DEFLATE), and overviews
- EPSG projection is read from HDF5 metadata (typically EPSG:3413 for polar scenes)
- Covariance terms include: HHHH, VVVV, HVHV, VHVH, HHVV, VVHH
- PySTAC items include raster metadata (units, datatype, EPSG, backscatter type)

## Running Tests

- Recommended: use conda-forge packages for GDAL/raster stack.
- Create env (example):
  - `conda create -n nisarhdf-dev -c conda-forge python=3.10 gdal rioxarray rasterio geopandas numpy pytest pip`
  - `conda activate nisarhdf-dev`
- Install this package in editable mode: `pip install -e .`
- Run integration tests that exercise VRT writing/reading:
  - Only VRT-on-TIFF: `pytest -q test/test_vrt_write_real_gdal.py::test_writeMultiBandVrt_tiff_inputs`
  - Full integration set: `pytest -q -m integration`
- GCOV test input:
  - Provide a real GCOV path via `--gcov /path/to/GCOV.h5` or env `NISAR_GCOV_PATH`.
  - Remote paths are supported if accessible (e.g., `https://...` or `s3://...`).

Notes
- GDAL 3.9+ may ignore `relativeToVRT` in some contexts; the code handles this by trying explicit options and falling back safely.
- Installing `gdal` via pip is not recommended; use conda-forge for consistent binaries.
