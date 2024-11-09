#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 10:44:13 2024

@author: ian
"""
from nisarhdf.nisarBaseRangeDopplerHDF import nisarBaseRangeDopplerHDF
import os
import numpy as np

from nisarhdf import writeMultiBandVrt, formatGeojson

class nisarRUNWHDF(nisarBaseRangeDopplerHDF):
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
        nisarBaseRangeDopplerHDF.__init__(self,
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
        self.lookType = 'ML'
        self.productParams = ['NumberRangeLooks', 'NumberAzimuthLooks']
        for param in self.RDParams:
            self.productParams.append(f'{self.lookType}{param}')
            

   
    def getSlantRangeData(self):
        ''' Get slant range data '''
        bands = self.h5[self.product][self.bands]
        self.slantRangeData = np.array(
            bands[self.frequency][self.productType]['slantRange'])

    def getZeroDopplerTimeData(self):
        ''' Get zero Doppler time '''
        bands = self.h5[self.product][self.bands]
        self.zeroDopplerTimeData = np.array(
            bands[self.frequency][self.productType]['zeroDopplerTime'])
            
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
        self.getSLCSlantRange()
        self.getSLCZeroDopplerTime()
        
        self.getMLSize()
        self.getMLSlantRange()
        
        self.effectivePRF()
        self.getRangeErrorCorrection()
        if not secondary:
            self.getMLZeroDopplerTime()
        else:
            self.ImageName = 'secondary'
            self.getMLZeroDopplerTime(secondary=True)
            #
        self.getCorrectedTime()
        # Note if secondary, it will have been passed the reference (see below)
        self.orbit = self.parseStateVectors(XMLOrbit=self.referenceOrbitXML)
        #
        self.getOrbitPassDirection()
        self.getCorners()
        self.getRangeErrorCorrection()
        self.getTimeToFirstSample()
        self.getSkewOffsets()
        self.getCenterIncidenceAngle()
        self.getSquint()
        self.getDeltaT()
        self.getCenterLatLon()
        self.getSceneCenterSatelliteHeight()
        #self.getOrbitPassDirection()
        self.getWavelength()
        #self.getSingleLookPixelSize()
        self.getExtent()
        #
        # If not secondary, create the secondary and parse
        if not secondary:
            self.secondary = \
                nisarRUNWHDF(isSecondary=True,
                             referenceOrbitXML=self.secondaryOrbitXML,
                             debug=self.debug)
            self.secondary.h5 = self.h5
            self.secondary.parseParams(secondary=True)
            #self.secondaryDatetime = self.secondary.datetime
            #self.secondaryDate = self.secondary.datetime
            self.dT = np.round((self.secondaryDatetime -
                               self.datetime).total_seconds()/86400.,
                               decimals=3)
            self.secondary.referenceOrbit = self.secondaryOrbit
            self.secondary.frame = self.frame
        self.genGeodatProperties()
        fields = ['coherenceMagnitude', 'connectedComponents',
                  'ionospherePhaseScreen', 'ionospherePhaseScreenUncertainty',
                  'unwrappedPhase']
        self.loadData(fields)

    # def SCLSceneCenterAlongTrackSpacing(self):
    #     '''
    #     Get SLC Scene Center Spacing

    #     Returns
    #     -------
    #     None.

    #     '''
    #     productType = \
    #         self.h5[self.product][self.bands][self.frequency][self.productType]
    #     MLAzimuthSize = np.array(
    #         productType['sceneCenterAlongTrackSpacing']).item()
    #     self.SLCAzimuthPixelSize = MLAzimuthSize / self.NumberRangeLooks

    # def writeData(self, filename, productField, includeVRT=True,
    #               secondary=True,
    #               includeGeojson=True,
    #               geojsonName=None, geojsonNameSecondary=None,
    #               dataType='>f4',
    #               metaData=None):
    #     '''

    #     Parameters
    #     ----------
    #     filename : str
    #         Name of file to write data to.
    #     productField : str
    #         Key for data (e.g., unwrapped phase)
    #     includeVRT : Bool, optional
    #         Include filename.vrt file. The default is True.
    #     includeGeodat : TYPE, optional
    #         Include GrIMP filename.geodat file. The default is True.
    #     dataType : str, optional
    #         Data type to save as. The default is '>f4'.

    #     Returns
    #     -------
    #     None.

    #     '''
    #     if not hasattr(self, productField):
    #         self.getImageData(productField)
    #     #
    #     # Write the image data to a floating point file
    #     self._writeImageData(filename, getattr(self, productField), dataType)
    #     #
    #     # write geodat in same dir as filename with geodatNRxNA unless other
    #     if includeGeojson:
    #         self.writeGeodatGeojson(filename=geojsonName,
    #                                 path=os.path.dirname(filename),
    #                                 secondary=secondary)
    #     #
    #     # Write an accompanyting vrt file if requested filename.vrt
    #     if includeVRT:
    #         if '>' in dataType:
    #             byteOrder = 'MSB'
    #         else:
    #             byteOrder = 'LSB'
    #         # pixel centered, lower left corner, pixel units
    #         geoTransform = [-0.5, 1., 0., -0.5, 0., 1.]
    #         writeMultiBandVrt(f'{filename}.vrt',
    #                           self.MLRangeSize, self.MLAzimuthSize,
    #                           [os.path.basename(filename)],
    #                           [productField],
    #                           metaData=metaData,
    #                           byteOrder=byteOrder,
    #                           geoTransform=geoTransform)
