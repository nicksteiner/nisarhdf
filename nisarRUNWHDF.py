#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:44:13 2024

@author: ian
"""
from nisarHDF import nisarBaseHDF
import os
import numpy as np

class nisarRUNWHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', product='RUNW', frequency='frequencyA',
                 productType='interferogram', polarization='HH', layer=None,
                 productData='unwrappedPhase', bands='swaths', isSecondary=False,
                 referenceOrbitXML=None, secondaryOrbitXML=None, debug=False):
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
                              bands=bands,
                              isSecondary=isSecondary,
                              referenceOrbitXML=referenceOrbitXML, 
                              secondaryOrbitXML=secondaryOrbitXML,
                              debug=debug)

    def parseParams(self):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        self.getGranuleNames()
        self.getOrbitAndFrame()
        self.getLookDirection()
        self.parseRefDate()
        self.getSlantRange()
        self.getZeroDopplerTime()
        self.getNumberOfLooks()
        self.getSize()
        self.effectivePRF()
        self.getRangeErrorCorrection()
        self.orbit = self.parseStateVectors(XMLOrbit=self.referenceOrbitXML)
        self.getCorners()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        self.getCenterIncidenceAngle()
        self.getSquint()
        self.getDeltaT()
        self.getCorrectedTime()
        self.getCenterLatLon()
        self.getSatelliteHeight()
        self.getOrbitPassDirection()
        self.getWavelength()
        self.getSingleLookPixelSize()
        #
        print('xxxxxxxx')
        self.secondary = nisarRUNWHDF(isSecondary=True, secondaryOrbitXML=self.secondaryOrbitXML, debug=self.debug)
        self.secondary.h5 = self.h5
        self.secondary.parseSecondary()

    def parseSecondary(self, secondaryOrbitXML=None):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        print('fffff')
        self.ImageName = 'secondary'
        self.getLookDirection()
        self.parseRefDate()
        self.getSlantRange()
        self.effectivePRF()
        self.getNumberOfLooks()
        self.getZeroDopplerTimeSecondary()
        self.getSize()
        self.getRangeErrorCorrection()
        print(self.secondaryOrbitXML)
        self.orbit = self.parseStateVectors(XMLOrbit=self.secondaryOrbitXML)
        self.getCorners()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        #
        self.getCenterIncidenceAngle()
        self.getSquint()
        self.getDeltaT()
        self.getCorrectedTime()
        self.getCenterLatLon()
        self.getSatelliteHeight()
        self.getOrbitPassDirection()
        self.getWavelength()
        self.getSingleLookPixelSize()

    def maskPhase(self, largest=True):
        '''
        Mask phase using connected components mask
        largest : Bool
            Retain only the largest connected component.
        Returns
        -------
        None.

        '''
        #
        # Ensure data have been extracted
        if not hasattr(self, 'connectedComponents'):
            self.getImageData('connectedComponents')
        if not hasattr(self, 'unwrappedPhase'):
            self.getImageData('unwrappedPhase')
        #
        #  Make a copy
        self.maskedUnwrappedPhase = self.unwrappedPhase.copy()
        # Eliminate phases with connected component values of 0
        self.maskedUnwrappedPhase[self.connectedComponents < 1] = np.nan
        #
        # Retain only the larged connected compnent.
        if largest:
            # Get the nonzero compponent labels
            labels = np.unique(self.connectedComponents)
            labels = labels[labels > 0]
            # Count for each label value
            count = np.array([np.sum(self.connectedComponents == x)
                              for x in labels])
            # find the max
            smallLabels = labels[count < np.max(count)]
            print('smallLabels', smallLabels)
            #
            # null out all the smaller
            for label in smallLabels:
                self.maskedUnwrappedPhase[label ==
                                          self.connectedComponents] = np.nan

    def writeData(self, filename, productField, includeVRT=True, 
                  secondary=True,
                  includeGeojson=True,
                  geojsonName=None, geojsonNameSecondary=None,
                  dataType='>f4',
                  metaData=None):
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
        # Write the image data to a floating point file
        self._writeImageData(filename, getattr(self, productField), dataType)
        #
        # write geodat in same dir as filename with geodatNRxNA unless other
        if includeGeojson:
            self.writeGeodatGeojson(filename=geojsonName,
                                    path=os.path.dirname(filename),
                                    secondary=secondary)
        #
        # Write an accompanyting vrt file if requested filename.vrt
        if includeVRT:
            if '>' in dataType:
                byteOrder = 'MSB'
            else:
                byteOrder = 'LSB'
            # pixel centered, lower left corner, pixel units
            geoTransform = [0.5, 0.5, 1., 0., 0., 1.]
            self._writeVrt(f'{filename}.vrt', [filename], [productField],
                           metaData=metaData, byteOrder=byteOrder,
                           setSRS=False, geoTransform=geoTransform)