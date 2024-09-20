#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:44:13 2024

@author: ian
"""
from nisarhdf import nisarBaseHDF
import os
import numpy as np
import rasterio


class nisarRUNWHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', frequency='frequencyA',
                 polarization='HH', isSecondary=False,
                 referenceOrbitXML=None, secondaryOrbitXML=None, debug=False):
        '''
       sar : str, optional
            SAR Idendtifier (always LSAR for now). The default is 'LSAR'.
        frequency : str, optional
            frequency band to extract. The default is 'frequencyA'.
        polarization : str, optional
            Polarization. The default is 'HH'.
        isSecondary : Boolean, optional
            For internal use only. The default is False.
        referenceOrbitXML : str, optional
            XML file to override orbit in hdf The default is None.
        secondaryOrbitXML : TYPE, optional
            XML file for secondary orbit. The default is None

        Returns
        -------
        None.

        '''
        print('ref orbit', referenceOrbitXML)
        nisarBaseHDF.__init__(self,
                              sar=sar,
                              product='RUNW',
                              frequency=frequency,
                              productType='interferogram',
                              polarization=polarization,
                              layer=None,
                              productData='unwrappedPhase',
                              bands='swaths',
                              isSecondary=isSecondary,
                              referenceOrbitXML=referenceOrbitXML,
                              secondaryOrbitXML=secondaryOrbitXML,
                              debug=debug)

    def parseParams(self, secondary=False, **keywords):
        '''
        Parse all the params needed to make a geodatNRxNA.geojson file

        Returns
        -------
        None.

        '''
        self.getGranuleNames()
        self.getOrbitAndFrame(**keywords)
        self.getLookDirection()
        self.parseRefDate()
        self.getNumberOfLooks()
        self.getSize()
        self.getSlantRange()

        self.effectivePRF()
        self.getRangeErrorCorrection()
        if not secondary:
            self.getZeroDopplerTime()
        else:
            self.ImageName = 'secondary'
            self.getZeroDopplerTimeSecondary()
            #
        self.getCorrectedTime()
        # Note if secondary, it will have been passed the reference (see below)
        self.orbit = self.parseStateVectors(XMLOrbit=self.referenceOrbitXML)
        #
        self.getCorners()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        self.getCenterIncidenceAngle()
        self.getSquint()
        self.getDeltaT()
        self.getCenterLatLon()
        self.getSatelliteHeight()
        self.getOrbitPassDirection()
        self.getWavelength()
        self.getSingleLookPixelSize()
        #
        # If not secondary, create the secondary and parse
        if not secondary:
            self.secondary = \
                nisarRUNWHDF(isSecondary=True,
                             referenceOrbitXML=self.secondaryOrbitXML,
                             debug=self.debug)
            self.secondary.h5 = self.h5
            self.secondary.parseParams(secondary=True)
            self.secondaryDatetime = self.secondary.datetime
            self.secondaryDate = self.secondary.datetime
            self.dT = np.round((self.secondaryDatetime -
                               self.datetime).total_seconds()/86400.,
                               decimals=3)

    def applyMask(self, maskFile):
        '''
        Apply an external mask (e.g., icemask) file

        Parameters
        ----------
        maskFile : str
            Path to mask file.

        Returns
        -------
        None.

        '''
        if not os.path.exists(maskFile):
            self.printError('applyMask: file ({maskFile} does not exist')
            return

        if not hasattr(self, 'maskedUnwrappedPhase'):
            self.maskedUnwrappedPhase = self.unwrappedPhase.copy()
        #
        mask = rasterio.open(maskFile).read()[0]
        self.maskedUnwrappedPhase[mask < 1] = np.nan

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
        if not hasattr(self, 'maskedUnwrappedPhase'):
            self.maskedUnwrappedPhase = self.unwrappedPhase.copy()
        # Eliminate phases with connected component values of 0
        print(np.sum(self.connectedComponents < 1))
        #
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
            geoTransform = [-0.5, 1., 0., -0.5, 0., 1.]
            self._writeVrt(f'{filename}.vrt',
                           [os.path.basename(filename)],
                           [productField],
                           metaData=metaData,
                           byteOrder=byteOrder,
                           setSRS=False,
                           geoTransform=geoTransform)
