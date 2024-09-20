#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:21:35 2024

@author: ian
"""

from nisarHDF import nisarBaseHDF
import numpy as np


class nisarGUNWHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', product='GUNW', frequency='frequencyA',
                 productType='unwrappedInterferogram', polarization='HH',
                 layer=None, productData='unwrappedPhase', bands='grids'):
        '''
        Instantiate nisarVel object. Possible bands are 'image', 'sigma0',
        'gamma0', or user defined.
        Parameters
        ----------
        imageType: str
            imageType custom name, or image, sigma0, gamma0. If not specified
            determined from grimp product name.
        verbose : bool, optional
            Increase level of informational messages. The default is True.
        noData : scalar
            no data value. Defaults to np.nan if not image/sigma/gamma
        Returns
        -------
        None.
        '''
        nisarBaseHDF.__init__(self,
                              sar=sar,
                              product=product,
                              frequency=frequency,
                              productType=productType,
                              polarization=polarization,
                              layer=layer,
                              productData=productData,
                              bands=bands)

    def parseParams(self):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        self.getOrbitAndFrame()
        self.getLookDirection()
        self.parseRefDate()
        self.getNumberOfLooks()
        self.getSize()
        self.getWavelength()
        self.getGeoCoordinates()

    def writeData(self, filename, productField, includeVRT=True,
                  includeGeodat=True, dataType='>f4', metaData=None):
        '''

        Parameters
        ----------
        filename : str
            Name of file to write data to.
        productField : str
            Key for data (e.g., unwrapped phase)
        includeVRT : Bool, optional
            Include filename.vrt file. The default is True.
        includeGeodat : TYPE, optional
            Include GrIMP filename.geodat file. The default is True.
        dataType : str, optional
            Data type to save as. The default is '>f4'.

        Returns
        -------
        None.

        '''
        if not hasattr(self, productField):
            self.getImageData(productField)
        #
        if self.dy < 0:
            data = np.flipud(getattr(self, productField))
        else:
            data = self.getattr(self, productField)
        #
        self._writeImageData(filename, data, dataType)
        #
        if includeGeodat:
            self.writeGeodat(f'{filename}.geodat')
        if includeVRT:
            if '>' in dataType:
                byteOrder = 'MSB'
            else:
                byteOrder = 'LSB'
            geoTransform = [self.x0, np.abs(self.dx), 0,
                            self.y0, 0, np.abs(self.dy)]
            self._writeVrt(f'{filename}.vrt', [filename], [productField],
                           metaData=metaData, byteOrder=byteOrder,
                           setSRS=True, geoTransform=geoTransform)
