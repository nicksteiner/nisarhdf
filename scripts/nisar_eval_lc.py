#!/usr/bin/env python3
import argparse
import json
import math
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform
from rasterio.windows import Window
from affine import Affine

try:
    import pandas as pd
except ImportError:
    pd = None

try:
    import geopandas as gpd
    from shapely.geometry import Point
except ImportError:
    gpd = None
    Point = None


WC21 = {
    10: ("Tree cover", "#006400"),
    20: ("Shrubland", "#ffbb22"),
    30: ("Grassland", "#ffff4c"),
    40: ("Cropland", "#f096ff"),
    50: ("Built-up", "#fa0000"),
    60: ("Bare / sparse vegetation", "#b4b4b4"),
    70: ("Snow and ice", "#f0f0f0"),
    80: ("Permanent water bodies", "#0064c8"),
    90: ("Herbaceous wetland", "#0096a0"),
    95: ("Mangroves", "#00cf75"),
    100: ("Moss and lichen", "#fae6a0"),
}


@dataclass
class Options:
    sar_path: str
    wc_path: str
    out_csv: str
    out_gpkg: Optional[str]
    samples_per_class: int
    classes: Optional[List[int]]
    sar_nodata: Optional[float]
    wc_nodata: Optional[int]
    random_seed: int
    ptiles: Tuple[float, float]


def parse_args() -> Options:
    ap = argparse.ArgumentParser(
        description="Sample SAR dB values by ESA WorldCover 2021 classes."
    )
    ap.add_argument("--sar", required=True, help="Path to SAR raster (dB, GeoTIFF)")
    ap.add_argument("--worldcover", required=True, help="Path to ESA WorldCover 2021 raster (GeoTIFF)")
    ap.add_argument("--out-csv", required=True, help="Output CSV for per-class stats")
    ap.add_argument("--out-gpkg", help="Optional GeoPackage for sampled points")
    ap.add_argument("--samples-per-class", type=int, default=5000, help="Max random samples per class (default 5000)")
    ap.add_argument("--classes", nargs="+", type=int, help="Subset of WC classes to sample (e.g. 10 20 50)")
    ap.add_argument("--sar-nodata", type=float, default=None, help="Override SAR NoData (if None, taken from raster)")
    ap.add_argument("--wc-nodata", type=int, default=None, help="Override WC NoData (if None, taken from raster)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--percentiles", nargs=2, type=float, default=(10.0, 90.0), help="Two percentiles, e.g. 10 90")
    args = ap.parse_args()
    return Options(
        sar_path=args.sar,
        wc_path=args.worldcover,
        out_csv=args.out_csv,
        out_gpkg=args.out_gpkg,
        samples_per_class=args.samples_per_class,
        classes=args.classes,
        sar_nodata=args.sar_nodata,
        wc_nodata=args.wc_nodata,
        random_seed=args.seed,
        ptiles=(args.percentiles[0], args.percentiles[1]),
    )


def reproject_worldcover_to_sar(wc_ds, sar_ds) -> np.ndarray:
    """Reproject WC raster to match SAR grid using nearest neighbor."""
    dst = np.empty((sar_ds.height, sar_ds.width), dtype=np.int32)
    dst.fill(-2147483648)  # a sentinel we'll treat as nodata for wc reprojection
    reproject(
        source=rasterio.band(wc_ds, 1),
        destination=dst,
        src_transform=wc_ds.transform,
        src_crs=wc_ds.crs,
        dst_transform=sar_ds.transform,
        dst_crs=sar_ds.crs,
        dst_nodata=-2147483648,
        resampling=Resampling.nearest,
        num_threads=2,
    )
    return dst


def random_choice_safe(idx_flat: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    n = idx_flat.size
    if n == 0:
        return idx_flat
    if k >= n:
        return idx_flat
    return rng.choice(idx_flat, size=k, replace=False)


def indices_to_xy(rows: np.ndarray, cols: np.ndarray, transform: Affine) -> Tuple[np.ndarray, np.ndarray]:
    xs, ys = rasterio.transform.xy(transform, rows, cols)
    return np.asarray(xs), np.asarray(ys)


def main():
    opts = parse_args()
    rng = np.random.default_rng(opts.random_seed)

    if pd is None:
        raise RuntimeError("pandas is required. pip install pandas")
    if opts.out_gpkg and gpd is None:
        raise RuntimeError("geopandas is required for --out-gpkg. pip install geopandas shapely pyproj fiona")

    # Open datasets
    with rasterio.open(opts.sar_path) as sar_ds, rasterio.open(opts.wc_path) as wc_ds:
        sar_nodata = sar_ds.nodata if opts.sar_nodata is None else opts.sar_nodata
        wc_nodata = wc_ds.nodata if opts.wc_nodata is None else opts.wc_nodata

        # Reproject WC to SAR grid if needed
        same_grid = (
            sar_ds.crs == wc_ds.crs and
            sar_ds.transform.almost_equals(wc_ds.transform) and
            sar_ds.width == wc_ds.width and sar_ds.height == wc_ds.height
        )

        if same_grid:
            wc = wc_ds.read(1)
        else:
            wc = reproject_worldcover_to_sar(wc_ds, sar_ds)

        sar = sar_ds.read(1).astype("float32")
        # Build valid mask
        valid = np.ones(sar.shape, dtype=bool)
        if sar_nodata is not None and not math.isnan(sar_nodata):
            valid &= sar != sar_nodata
        valid &= np.isfinite(sar)

        if wc_nodata is not None:
            valid &= wc != wc_nodata
        # also drop sentinel from reprojection if any
        valid &= wc != -2147483648

        # choose classes
        all_classes = np.array(sorted(list(WC21.keys())), dtype=np.int32)
        if opts.classes:
            classes = np.array(sorted([c for c in opts.classes if c in WC21]), dtype=np.int32)
        else:
            # only keep classes that exist in the raster to speed things up
            present = np.unique(wc[valid])
            classes = np.array([c for c in all_classes if c in present], dtype=np.int32)

        rows_all, cols_all = np.nonzero(valid)
        # Mapping (row,col) -> flat index for faster subsetting per class
        flat_idx = rows_all * sar.shape[1] + cols_all
        wc_flat = wc[rows_all, cols_all]
        sar_flat = sar[rows_all, cols_all]

        # gather results
        records = []
        sample_points = []  # (x, y, class_code, class_name, sar_db)

        for cls in classes:
            # mask class
            mask_cls = wc_flat == cls
            idx_cls = flat_idx[mask_cls]
            if idx_cls.size == 0:
                continue

            pick = random_choice_safe(idx_cls, opts.samples_per_class, rng)
            # back to rows/cols
            rows = pick // sar.shape[1]
            cols = pick % sar.shape[1]
            vals = sar[rows, cols]

            # robust stats
            if vals.size == 0:
                continue
            p_lo, p_hi = np.percentile(vals, [opts.ptiles[0], opts.ptiles[1]])
            rec = {
                "class_code": int(cls),
                "class_name": WC21[int(cls)][0],
                "n": int(vals.size),
                "mean_db": float(np.mean(vals)),
                "median_db": float(np.median(vals)),
                "std_db": float(np.std(vals, ddof=1)) if vals.size > 1 else float("nan"),
                f"p{int(opts.ptiles[0])}_db": float(p_lo),
                f"p{int(opts.ptiles[1])}_db": float(p_hi),
            }
            records.append(rec)

            if opts.out_gpkg:
                xs, ys = indices_to_xy(rows, cols, sar_ds.transform)
                for x, y, v in zip(xs, ys, vals):
                    sample_points.append((x, y, int(cls), WC21[int(cls)][0], float(v)))

        # write CSV
        df = pd.DataFrame.from_records(records).sort_values("class_code")
        df.to_csv(opts.out_csv, index=False)
        print(f"Wrote stats: {opts.out_csv}")

        # write points (optional)
        if opts.out_gpkg and sample_points:
            gdf = gpd.GeoDataFrame(
                {
                    "class_code": [r[2] for r in sample_points],
                    "class_name": [r[3] for r in sample_points],
                    "sar_db": [r[4] for r in sample_points],
                },
                geometry=[Point(r[0], r[1]) for r in sample_points],
                crs=sar_ds.crs,
            )
            layername = os.path.splitext(os.path.basename(opts.out_gpkg))[0]
            gdf.to_file(opts.out_gpkg, layer=layername, driver="GPKG")
            print(f"Wrote samples: {opts.out_gpkg} ({len(gdf)} points)")
        elif opts.out_gpkg:
            print("No points written (no valid samples found).")


if __name__ == "__main__":
    main()
