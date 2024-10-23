#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 08:25:27 2024

@author: ian
"""
from osgeo import gdal, osr
import os


def writeMultiBandVrt(newVRTFile, xSize, ySize, sourceFiles, descriptions,
                      byteOrder=None, eType=gdal.GDT_Float32,
                      geoTransform=[-0.5, 1., 0., -0.5, 0., 1.], metaData=None,
                      epsg=None, noDataValue=-2.0e9):
    '''
    Write a vrt for for several products as bands. 
    Parameters
    ----------
    newVRTFile : str
        Name for vrt file being produced.
    sourceFiles : list
        List of files to serve as bands in vrt.
    descriptions : list
        List of descriptions for the bands.
    byteOrder : str, optional
        Byte order (MSB, LSB). If none, the default is use value from the
        metaData, if none there default to MSMB. The default is None.
    eType : GDT, optional
        Gdal data type. The default is gdal.GDT_Float32.
    geoTransform : list, optional
        The geoTransform. The default is [-0.5, 1., 0., -0.5, 0., 1.] for
        array with no coordinates.
    metaData : dict, optional
        Dictionary with meta data values. The default is None.
    epsg : str, optional
        Seet to epsg value to set SRS. The default is None.
    noDataValue : int, float as appropriate, optional
        The nodata value. The default is -2.0e9.

    Returns
    -------
    None.

    '''
    if type(sourceFiles) is not list:
        sourceFiles = [sourceFiles]
    if type(descriptions) is not list:
        descriptions = [descriptions]
    if metaData is None:
        metaData = {}
    # Kill any old file
    if os.path.exists(newVRTFile):
        os.remove(newVRTFile)
    # Create VRT
    bands = len(sourceFiles)
    drv = gdal.GetDriverByName("VRT")
    vrt = drv.Create(newVRTFile, xSize, ySize, bands=0, eType=eType)
    vrt.SetGeoTransform(geoTransform)
    #
    if epsg is not None:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        vrt.SetSpatialRef(sr)
    # Set byte order, note defaults to MSB for GrIMP
    if byteOrder is None:
        if "ByteOrder" in metaData:
            byteOrder = metaData["ByteOrder"]
        else:
            byteOrder = "MSB"
            metaData["ByteOrder"] = byteOrder
    else:
        metaData["ByteOrder"] = byteOrder
    # Set met data
    if metaData is not None:
        vrt.SetMetadata(metaData)
    # Loop to add bands
    for sourceFile, description, bandNumber in \
            zip(sourceFiles, descriptions, range(1, bands + 1)):
        options = [f"SourceFilename={os.path.basename(sourceFile)}",
                   "relativeToVRT=1",
                   "subclass=VRTRawRasterBand",
                   f"BYTEORDER={byteOrder}",
                   bytes(0)]
        vrt.AddBand(eType, options=options)
        band = vrt.GetRasterBand(bandNumber)
        band.SetMetadataItem("Description", description)
        band.SetNoDataValue(noDataValue)
    # Close the vrt
    vrt = None