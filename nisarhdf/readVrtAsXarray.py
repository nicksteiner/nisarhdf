#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:02:52 2024

@author: ian
"""
import rioxarray


def readVrtAsXarray(filename, mask_and_scale=True):
    '''
    Read a vrt and save each band as named variabile. Note this requires the
    band name be save and metadata field "Description" for each band as is
    for vrts saved with nisarHDF functions

    Parameters
    ----------
    filename : str
        Vrt file.

    Returns
    -------
    xarray : rioxarray.
        Data as a rioxarray.
    '''
    xarray = rioxarray.open_rasterio(
        filename,
        band_as_variable=True,
        mask_and_scale=mask_and_scale,
    )
    # Prefer GDAL band Description (and common alternatives); fallback to existing names
    rename_map = {}
    preferred_keys = (
        "Description",
        "description",
        "DESCRIPTION",
        "long_name",
        "standard_name",
        "name",
    )
    for var_name in list(xarray.data_vars):
        attrs = getattr(xarray[var_name], "attrs", {}) or {}
        new_name = None
        for k in preferred_keys:
            val = attrs.get(k)
            if isinstance(val, str) and len(val) > 0:
                new_name = val
                break
        if new_name is None:
            new_name = var_name
        rename_map[var_name] = new_name
    if any(rename_map[k] != k for k in rename_map):
        xarray = xarray.rename(rename_map)
    return xarray

