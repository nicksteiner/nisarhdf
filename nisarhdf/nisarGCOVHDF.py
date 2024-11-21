#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 08:51:41 2024

@author: ian
"""

from nisarhdf import nisarBaseGeocodedHDF


class nisarGCOVHDF(nisarBaseGeocodedHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', product='GCOV', frequency='frequencyA',
                 productType=None, polarization=None,
                 layer=None, productData=None, bands='grids',
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
        self.productParams = ['NumberRangeLooks', 'NumberAzimuthLooks']
        self.lookType = None

    def parseParams(self, secondary=False, noLoadData=False, **keywords):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        print(self.layer)
        self.getOrbitAndFrame(**keywords)
        self.getNumberOfLooks()
        self.getLookDirection()
        self.getOrbitPassDirection()
        self.parseRefDate()
        self.getGeoCoordinates()
        self.getWavelength()
        self.getGranuleNames()
        self.getEPSG()

        self.getSLCSlantRange()
        self.getSLCZeroDopplerTime()
        self.effectivePRF()
        self.getExtent()
        fields = [self.parseString(x) for x in
                  self.h5[self.product][self.bands][self.frequency][
                      'listOfCovarianceTerms']]
        #
        fields += ['mask', 'numberOfLooks', 'rtcGammaToSigmaFactor']
        self.loadData(fields, noLoadData=noLoadData)
