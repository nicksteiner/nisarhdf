#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:21:35 2024

@author: ian
"""

from nisarhdf import nisarBaseGeocodedHDF
import numpy as np
from nisarhdf import writeMultiBandVrt

class nisarGOFFHDF(nisarBaseGeocodedHDF):
    '''
    This class creates objects to work with nisar ROFFHDF images.
    '''

    def __init__(self,  sar='LSAR', product='GOFF', frequency='frequencyA',
                 productType='pixelOffsets', polarization='HH',
                 layer='layer1', productData='alongTrackOffset', bands='grids',
                 referenceOrbitXML=None, secondaryOrbitXML=None,
                 isSecondary=False):
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
        nisarBaseGeocodedHDF.__init__(self,
                                      sar=sar,
                                      product=product,
                                      frequency=frequency,
                                      productType=productType,
                                      polarization=polarization,
                                      layer=layer,
                                      productData=productData,
                                      bands=bands,
                                      isSecondary=isSecondary,
                                      referenceOrbitXML=referenceOrbitXML,
                                      secondaryOrbitXML=secondaryOrbitXML)
        self.productParams = ['xSize', 'ySize', 'dx', 'dy']

    def parseParams(self, secondary=False, **keywords):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        self.getOrbitAndFrame(**keywords)
        self.getLookDirection()
        self.getOrbitPassDirection()
        self.parseRefDate()
        self.getGeoCoordinates()
        self.getWavelength()
        self.effectivePRF()
        self.getOffsetWindowParams()
        self.getGranuleNames()
        self.getEPSG()
        self.getSLCSlantRange()
        self.getSLCZeroDopplerTime()
        self.getExtent()
        self.fieldDict = {'pixelOffsets': ['slantRangeOffset',
                                           'slantRangeOffsetVariance',
                                           'alongTrackOffset',
                                           'alongTrackOffsetVariance',
                                           'crossOffsetVariance',
                                           'correlationSurfacePeak',
                                           'snr']}
        self.getLayers(self.fieldDict['pixelOffsets'])
        self.scaleFactors = {'slantRangeOffset': 1./self.SLCRangePixelSize,
                             'alongTrackOffset': 1./self.SLCAzimuthPixelSize,
                             'slantRangeOffsetVariance':
                                 (1./self.SLCRangePixelSize)**2,
                             'alongTrackOffsetVariance':
                                 (1./self.SLCAzimuthPixelSize)**2}
        self.suffixes = {'slantRangeOffset': '.dr',
                         'slantRangeOffsetVariance': '.vr',
                         'alongTrackOffset': '.da',
                         'alongTrackOffsetVariance': '.va',
                         'crossOffsetVariance': '.vc',
                         'correlationSurfacePeak': '.cc',
                         'snr': '.snr'}
        