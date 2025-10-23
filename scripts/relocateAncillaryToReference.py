#!/usr/bin/env python3
"""
Relocate ancillary raster files to match the extent and spacing of a reference TIFF file.

This script resamples and reprojects ancillary data (DEM, landcover, etc.) to align
with a reference raster file's projection, extent, and pixel spacing.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Optional, Union

import numpy as np
from osgeo import gdal, osr

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('relocate_ancillary')


def get_raster_metadata(raster_path: str) -> dict:
    """
    Extract metadata from a raster file.
    
    Parameters
    ----------
    raster_path : str
        Path to the raster file.
    
    Returns
    -------
    dict
        Dictionary containing geotransform, projection, width, and length.
    """
    ds = gdal.Open(raster_path, gdal.GA_ReadOnly)
    if ds is None:
        raise FileNotFoundError(f"Unable to open raster file: {raster_path}")
    
    geotransform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    width = ds.RasterXSize
    length = ds.RasterYSize
    
    ds = None
    
    return {
        'geotransform': geotransform,
        'projection': projection,
        'width': width,
        'length': length
    }


def relocate_raster(
    input_path: str,
    output_path: str,
    reference_metadata: dict,
    resampling_method: str = 'near',
    no_data: Optional[float] = None,
    cog: bool = True
) -> None:
    """
    Relocate (resample and reproject) a raster to match reference metadata.
    
    Parameters
    ----------
    input_path : str
        Path to the input ancillary raster file.
    output_path : str
        Path to the output relocated raster file.
    reference_metadata : dict
        Metadata from the reference raster (geotransform, projection, width, length).
    resampling_method : str, optional
        Resampling method. Options: 'near', 'bilinear', 'cubic', 'cubicspline',
        'lanczos', 'average', 'mode'. Default is 'near'.
    no_data : float, optional
        No-data value to use for the output raster.
    cog : bool, optional
        Whether to create a Cloud Optimized GeoTIFF. Default is True.
    """
    # Map resampling method names to GDAL constants
    resampling_dict = {
        'near': gdal.GRA_NearestNeighbour,
        'bilinear': gdal.GRA_Bilinear,
        'cubic': gdal.GRA_Cubic,
        'cubicspline': gdal.GRA_CubicSpline,
        'lanczos': gdal.GRA_Lanczos,
        'average': gdal.GRA_Average,
        'mode': gdal.GRA_Mode
    }
    
    if resampling_method not in resampling_dict:
        raise ValueError(f"Invalid resampling method: {resampling_method}. "
                        f"Choose from {list(resampling_dict.keys())}")
    
    resampling_alg = resampling_dict[resampling_method]
    
    # Open source dataset
    src_ds = gdal.Open(input_path, gdal.GA_ReadOnly)
    if src_ds is None:
        raise FileNotFoundError(f"Unable to open input file: {input_path}")
    
    # Get source data type
    src_band = src_ds.GetRasterBand(1)
    src_dtype = src_band.DataType
    
    # Determine output data type and creation options
    if cog:
        creation_options = [
            'TILED=YES',
            'COMPRESS=LZW',
            'COPY_SRC_OVERVIEWS=YES',
            'GEOTIFF_VERSION=1.1'
        ]
    else:
        creation_options = ['COMPRESS=LZW']
    
    # Create output dataset
    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(
        output_path,
        reference_metadata['width'],
        reference_metadata['length'],
        src_ds.RasterCount,
        src_dtype,
        options=creation_options
    )
    
    if dst_ds is None:
        src_ds = None
        raise RuntimeError(f"Unable to create output file: {output_path}")
    
    # Set output geotransform and projection
    dst_ds.SetGeoTransform(reference_metadata['geotransform'])
    dst_ds.SetProjection(reference_metadata['projection'])
    
    # Set no-data value if provided
    if no_data is not None:
        for band_idx in range(1, dst_ds.RasterCount + 1):
            dst_band = dst_ds.GetRasterBand(band_idx)
            dst_band.SetNoDataValue(float(no_data))
    
    # Perform the reprojection/resampling
    logger.info(f"Reprojecting {input_path} to match reference...")
    logger.info(f"  Resampling method: {resampling_method}")
    logger.info(f"  Output size: {reference_metadata['width']} x {reference_metadata['length']}")
    
    gdal.ReprojectImage(
        src_ds,
        dst_ds,
        src_ds.GetProjection(),
        reference_metadata['projection'],
        resampling_alg
    )
    
    # Build overviews for COG if requested
    if cog:
        logger.info("Building overviews for COG...")
        dst_ds.BuildOverviews("NEAREST", [2, 4, 8, 16])
    
    # Close datasets
    src_ds = None
    dst_ds = None
    
    logger.info(f"Successfully created: {output_path}")


def main():
    """Main function to handle command-line interface."""
    parser = argparse.ArgumentParser(
        description='Relocate ancillary raster files to match a reference TIFF file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Relocate a DEM to match a reference image
  %(prog)s -r reference.tif -i dem.tif -o dem_relocated.tif
  
  # Relocate with bilinear resampling and custom no-data value
  %(prog)s -r reference.tif -i landcover.tif -o landcover_relocated.tif \\
           --method bilinear --nodata 255
  
  # Relocate without creating COG
  %(prog)s -r reference.tif -i input.tif -o output.tif --no-cog
        """
    )
    
    parser.add_argument(
        '-r', '--reference',
        required=True,
        help='Path to reference TIFF file (defines target extent and spacing)'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Path to input ancillary raster file to relocate'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output relocated raster file'
    )
    
    parser.add_argument(
        '-m', '--method',
        default='near',
        choices=['near', 'bilinear', 'cubic', 'cubicspline', 'lanczos', 'average', 'mode'],
        help='Resampling method (default: near)'
    )
    
    parser.add_argument(
        '--nodata',
        type=float,
        default=None,
        help='No-data value for output raster'
    )
    
    parser.add_argument(
        '--no-cog',
        action='store_true',
        help='Do not create Cloud Optimized GeoTIFF (COG)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Validate input files
    if not os.path.isfile(args.reference):
        logger.error(f"Reference file not found: {args.reference}")
        sys.exit(1)
    
    if not os.path.isfile(args.input):
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    try:
        # Get reference metadata
        logger.info(f"Reading reference metadata from: {args.reference}")
        ref_metadata = get_raster_metadata(args.reference)
        
        # Relocate the raster
        relocate_raster(
            input_path=args.input,
            output_path=args.output,
            reference_metadata=ref_metadata,
            resampling_method=args.method,
            no_data=args.nodata,
            cog=not args.no_cog
        )
        
        logger.info("Processing complete!")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
