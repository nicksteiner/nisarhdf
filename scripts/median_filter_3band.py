#!/usr/bin/env python3
import argparse
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.env import Env
from scipy.ndimage import median_filter, generic_filter

def parse_args():
    ap = argparse.ArgumentParser(
        description="Median filter a 3-band (or multiband) raster using SciPy."
    )
    ap.add_argument("inp", help="Input GeoTIFF")
    ap.add_argument("out", help="Output COG")
    ap.add_argument("--size", type=int, default=3, help="Kernel size (odd, e.g., 3,5)")
    ap.add_argument("--ignore-nodata", action="store_true",
                    help="Ignore NoData in window using nanmedian (slower but safer).")
    ap.add_argument("--threads", type=str, default="ALL_CPUS",
                    help="GDAL NUM_THREADS (default ALL_CPUS)")
    return ap.parse_args()

def nanmedian_func(win):
    # generic_filter callback: compute median ignoring NaNs
    return np.nanmedian(win)

def filter_band(arr, size, ignore_nodata, nodata):
    """Apply median filter on a single band array."""
    if ignore_nodata:
        # Convert nodata to NaN, use generic_filter with nanmedian
        work = arr.astype("float64", copy=False)
        if nodata is not None:
            work = work.copy()
            work[arr == nodata] = np.nan
        filt = generic_filter(work, nanmedian_func, size=(size, size), mode="nearest")
        # Restore dtype; keep nodata where center pixel was nodata
        if nodata is not None:
            center_mask = (arr == nodata)
            filt[center_mask] = nodata
        return filt.astype(arr.dtype, copy=False)
    else:
        # Fast path: standard median_filter; does not ignore nodata in window.
        # We preserve nodata at center so nodata areas don't get invented.
        center_mask = (arr == nodata) if nodata is not None else None
        filt = median_filter(arr, size=size, mode="nearest")
        if center_mask is not None:
            filt[center_mask] = nodata
        return filt

def main():
    args = parse_args()
    k = args.size if args.size % 2 == 1 else args.size + 1

    cog_profile = dict(
        driver="COG",
        compress="LZW",
        blockxsize=512,
        blockysize=512,
        BIGTIFF="YES",
        NUM_THREADS=args.threads,
        overview_resampling=Resampling.nearest.name,  # safe default
    )

    with Env(NUM_THREADS=args.threads):
        with rasterio.open(args.inp) as src:
            profile = src.profile
            nodata = src.nodata
            count = src.count
            # Read all bands (works for 3-band or more)
            data = src.read()  # shape: (bands, rows, cols)

        # Filter each band independently
        out = np.empty_like(data)
        for b in range(count):
            out[b] = filter_band(data[b], k, args.ignore_nodata, nodata)

        # Write COG
        profile.update(cog_profile)
        with rasterio.open(args.out, "w", **profile) as dst:
            dst.write(out)
            if nodata is not None:
                dst.update_tags(ns="rio_overview", resampling=Resampling.nearest.value)

    print(f"Median-filtered COG written: {args.out}")
    print(f"Bands processed: {count}, kernel: {k}x{k}, ignore_nodata={args.ignore_nodata}")

if __name__ == "__main__":
    main()
