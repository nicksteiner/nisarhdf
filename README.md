# nisarHDF

This repository contains python code to read the NISAR HDF formats and parse out the necessary parameters and reformatting necessary for the Greenland Ice Mapping Project (*GrIMP*) workflow. Most of the functionality though is independent of *GrIMP* and may prove useful for anyone working with NISAR HDF data.

In particular, it simplifies access to NISAR data by allowing:
  - Access data (e.g., `myRUNW.unwrappedPhase`) as numpy variables and other commonly used parameters (e.g., `myGUNW.NumberRangeLooks`) rather than hdf structures,
  - Mask data (e.g., remove  data with connected component values of zero),
  - Apply ionospheric corrections to phase data,
  - Interpolate metadata cubes (e.g., incidenceAngle and baseline), including mapping these variables to the same grid as the data products.
  - Save data and interpolated metadata cubes to binary or geotiff files, and
  - Write a .vrt file to link multiple data layers so they can be read back in as `xarray`.

At present classes have been developed for **RSLC**, **ROFF**, **RUNW**, **RIFG**, **GOFF**, **GUNW**, and **GCOV** products. A tutorial notebook for each of these data types is included in the *Notebook* folder in this repository.
