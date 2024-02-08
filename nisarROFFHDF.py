#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 08:36:38 2024

@author: ian
"""
from nisarHDF import nisarBaseHDF
import os
import numpy as np


class nisarROFFHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', product='ROFF', frequency='frequencyA',
                 productType='pixelOffsets', polarization='HH', layer='layer1',
                 productData='alongTrackOffset', bands='swaths'):
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
                              productData=productData)
    
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
        self.getSlantRange()
        self.getZeroDopplerTime()
        #self.getNumberOfLooks()
        self.getSize()
        self.effectivePRF(offsets=True)
        self.getRangeErrorCorrection()
        self.getCorners()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        self.getOffsetParams()
        self.getOffsetLayerParams()
        #self.getCenterIncidenceAngle()
        #self.getSquint()
        #self.getDeltaT()
        #self.getCorrectedTime()
        self.parseStateVectors()
        self.getCenterLatLon()
        self.getSatelliteHeight()
        self.getOrbitPassDirection()
        self.getWavelength()
        #self.getSingleLookPixelSize()