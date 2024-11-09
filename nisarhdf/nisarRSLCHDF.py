#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 08:05:24 2024

@author: ian
"""

from nisarhdf import nisarBaseHDF
from osgeo import gdal
import os
import numpy as np
import rioxarray
from nisarhdf import writeMultiBandVrt


class nisarRSLCHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', frequency='frequencyA', polarization='HH',
                 byteOrder='LSB'):
        '''
        NISAR ROFF HDF Reader

        Parameters
        ----------
        sar : str, optional
            SAR Idendtifier (always LSAR for now). The default is 'LSAR'.
        frequency : str, optional
            frequency band to extract. The default is 'frequencyA'.
        polarization : str, optional
            Polarization. The default is 'HH'.

        Returns
        -------
        None.

        '''
        nisarBaseHDF.__init__(self,
                              sar=sar,
                              product='RSLC',
                              frequency=frequency,
                              productType='pixelOffsets',
                              polarization=polarization,
                              # layer=layer,
                              productData=None,
                              bands='swaths',
                              bytOrder=byteOrder)
        self.isSecondary = False

    def parseParams(self, **kwds):
        '''
        Parse all the params needed for offsets

        Returns
        -------
        None.

        '''

        # self.getOffsetParams()
        # self.getNumberOfLooks()
        self.getDateTime()
        self.getSLCSLantRangeAndZeroDoppler()
        self.orbit = self.parseStateVectors(XMLOrbit=self.referenceOrbitXML,
                                            SLC=True)
        self.getSceneCenterSatelliteHeight()
        self.getSLCSize(**kwds)
        # self.getZeroDopplerTime()
        # self.getSize(offsets=True)
        # self.getSingleLookPixelSizeOffsets()
        # self.getEPSG()
        # self.effectivePRF()

    def writeSLC(self, fileName, frequency='frequencyA', dataType='>c8'):
        '''
        Parameters
        ----------
        fileName : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        swaths = self.h5Full['science'][self.sar][self.product]['swaths']
        data = np.array(swaths['frequencyA']['HH'].astype('c8'))
        self._writeImageData(fileName, data, dataType, noData=-2.e9)

    def writeOldParFile(self, filename):
        '''
        Write an old gamma style par file for SLC
        '''
        with open(filename, 'w') as fp:
            print(os.path.basename(self.hdfFile), file=fp)
            print(self.datetime.strftime('%d %b %Y').upper(), file=fp)
            print(self.datetime.strftime('%H %M %S.%f'), file=fp)
            print(' ', file=fp)
            print('0. 0.\n0.', file=fp)
            print(self.SpaceCraftAltitude, file=fp)
            print('0.0\n 0.0 0.0 0.0\n 0.0 0.0 0.0\n 0.0 0.0 0.0', file=fp)
            print(self.PRF, file=fp)
            print('0.0 0.0 0.0 0.0\n0.\n0.\n0. 0. 0.', file=fp)
            print(f'{self.SLCNearRange:0.2f}  {self.SLCCenterRange:.2f} '
                  f'{self.SLCFarRange:.2f}', file=fp)
            print(f'{self.SLCRangePixelSize:.8f}  0', file=fp)
            print('0\n1\n0\n0\n0\n1\n1 1\n0.', file=fp)
            print(f'{self.SLCAzimuthPixelSize:.8f} 0.', file=fp)
            print(f'{self.SLCRangeSize} {self.SLCAzimuthSize}', file=fp)
            print(f'{self.orbit.StateVectorInterval}', file=fp)
            print(f'{self.orbit.TimeOfFirstStateVector}', file=fp)
            
            for p, v in zip(self.orbit.position, self.orbit.velocity):
                print(f'{p[0]:10e} {p[1]:10e} {p[2]:10e}', file=fp)
                print(f'{v[0]:e} {v[1]:e} {v[2]:e}', file=fp)
                
                
        

    def getSLCSize(self, frequency='frequencyA', **kwds):
        '''
        Get size of SLC image

        Returns
        -------
        None.

        '''
        swath = self.h5Full['science'][self.sar][self.product]['swaths']
        self.SLCAzimuthSize, self.SLCRangeSize = \
            swath[frequency][self.polarization].shape

    def getDateTime(self):
        '''
        Get the time
        Returns
        -------
        None.

        '''
        id = self.h5Full['science'][self.sar]['identification']
        dateStr = np.array(id["zeroDopplerStartTime"]).item().decode('ascii')
        self.datetime, self.Date = self.parseDateStr(dateStr)

    def getSLCSLantRangeAndZeroDoppler(self, frequency='frequencyA'):
        '''
        Get near, center, and far slant range

        Returns
        -------
        None.

        '''
        swaths = self.h5Full['science'][self.sar][self.product]['swaths']
        self.SLCNearRange = swaths[frequency]['slantRange'][0]
        self.SLCFarRange = swaths[frequency]['slantRange'][-1]
        self.SLCCenterRange = (self.SLCNearRange + self.SLCFarRange) * 0.5
        #
        self.firstZeroDopplerTime = swaths['zeroDopplerTime'][0]
        self.lastZeroDopplerTime = swaths['zeroDopplerTime'][-1]
        self.midZeroDopplerTime = \
            (self.firstZeroDopplerTime + self.lastZeroDopplerTime) * 0.5
        self.PRF = 1. / np.array(swaths['zeroDopplerTimeSpacing']).item()
        #
        self.SLCRangePixelSize = np.array(
            swaths[frequency]['slantRangeSpacing']).item()
        self.SLCAzimuthPixelSize = np.array(
            swaths[frequency]['sceneCenterAlongTrackSpacing']).item()
