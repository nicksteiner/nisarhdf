#!/usr/bin/env python3
"""
Example Prefect flow wrapping nisar_gcov_cog.py.
Demonstrates structured logging integration and task-level observability.
"""
from pathlib import Path
from typing import List, Optional

from prefect import flow, task, get_run_logger


@task(name="convert-gcov-to-cog", retries=2)
def convert_gcov_to_cog(
    h5_path: Path,
    output_dir: Path,
    stac_dir: Path,
    sigma0: bool = False,
    db: bool = False,
    no_incidence_angle: bool = False,
    overwrite: bool = False,
    downsample_factor: int = 1,
) -> bool:
    """
    Task wrapper for nisar_gcov_cog main logic.
    Prefect's logger automatically captures all logging from the script.
    """
    logger = get_run_logger()
    logger.info("Converting GCOV %s to COG (sigma0=%s, db=%s)", h5_path.name, sigma0, db)

    # Import here to avoid module-level side effects; Prefect manages the logger tree
    from nisar_gcov_cog import main

    # Build argv-style args
    argv = [
        str(h5_path),
        "-o", str(output_dir),
        "--stac", str(stac_dir),
        "--downsample-factor", str(downsample_factor),
    ]
    if sigma0:
        argv.append("--sigma0")
    if db:
        argv.append("--db")
    if no_incidence_angle:
        argv.append("--no-incidence-angle")
    if overwrite:
        argv.append("--overwrite")

    try:
        exit_code = main(argv)
        if exit_code != 0:
            logger.error("Conversion failed with exit code %d", exit_code)
            return False
        logger.info("Conversion complete for %s", h5_path.name)
        return True
    except Exception as exc:
        logger.exception("Conversion raised exception for %s: %s", h5_path.name, exc)
        return False


@flow(name="gcov-to-cog-batch", log_prints=True)
def gcov_to_cog_batch(
    input_files: List[Path],
    output_dir: Path,
    stac_dir: Path,
    sigma0: bool = False,
    db: bool = False,
    no_incidence_angle: bool = False,
    overwrite: bool = False,
    downsample_factor: int = 1,
):
    """
    Batch convert multiple GCOV HDF5 files to COG with STAC metadata.
    Prefect tracks each file's task separately with full log capture.
    """
    logger = get_run_logger()
    logger.info("Starting batch conversion of %d files", len(input_files))

    output_dir.mkdir(parents=True, exist_ok=True)
    stac_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for h5_path in input_files:
        success = convert_gcov_to_cog(
            h5_path=h5_path,
            output_dir=output_dir,
            stac_dir=stac_dir,
            sigma0=sigma0,
            db=db,
            no_incidence_angle=no_incidence_angle,
            overwrite=overwrite,
            downsample_factor=downsample_factor,
        )
        results.append((h5_path.name, success))

    successful = sum(1 for _, s in results if s)
    logger.info("Batch complete: %d/%d succeeded", successful, len(results))

    return results


if __name__ == "__main__":
    # Example: process all .h5 in a directory
    input_dir = Path("/media/nsteiner/data3/nisar/oasis")
    output_dir = Path("./gcov_cogs")
    stac_dir = Path("./stac_items")

    h5_files = list(input_dir.glob("NISAR_L2_PR_GCOV_*.h5"))

    gcov_to_cog_batch(
        input_files=h5_files,
        output_dir=output_dir,
        stac_dir=stac_dir,
        db=True,
        no_incidence_angle=True,
        overwrite=False,
    )
