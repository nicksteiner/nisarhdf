#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 08:47:03 2024

@author: ian
"""
from abc import ABCMeta, abstractmethod
import h5py
import scipy
from datetime import datetime, timedelta
import numpy as np
import geojson
import pyproj
from osgeo import gdal, osr
from scipy import optimize
import os
from nisarHDF import nisarOrbit


class nisarBaseHDF():
    '''
    Abstract class to define parser for NISAR HDF products.
    '''

    __metaclass__ = ABCMeta

    def __init__(self, **keywords):
        '''
        Initialize a nisar hdf parser
        Parameters
        ----------
        sar : TYPE, optional
            hdf key. The default is 'LSAR'.
        product : TYPE, optional
            hdf key. The default is 'RUNW'.
        frequency : TYPE, optional
            hdf key.. The default is 'frequencyA'.
        productType : TYPE, optional
            hdf key.. The default is 'interferogram'.
        polarization : TYPE, optional
            hdf key.. The default is 'HH'.
        layer : TYPE, optional
            hdf key.. The default is None.
        productData : TYPE, optional
            Dhdf key.. The default is 'unwrappedPhase'.
        Returns
        -------
        None.

        '''
        print(keywords)
        self.initTags(**keywords)
        # Constants WGS84
        self.EarthRadiusMajor = 6378137.0
        self.EarthRadiusMinor = 6356752.3142
        self.EarthF = 1./298.257223563
        self.f = (self.EarthRadiusMajor -
                  self.EarthRadiusMinor) / self.EarthRadiusMajor
        self.cLight = 2.99792458e8
        self.epsg = None

    @abstractmethod
    def parseParams(**keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass

    @abstractmethod
    def writeData():
        ''' Abstract pointPlot method - define in child class. '''
        pass

    def initTags(self, **kwargs):
        '''
        Initialize tags as attributes

        Parameters
        ----------
        **kwargs : keywords
            Keywords to initialize or overwrite.

        Returns
        -------
        None.

        '''
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def printError(self, msg):
        '''
        Print error message
        Parameters
        ----------
        msg : str
            error message.
        Returns
        -------
        None

        '''
        print(f'\n\t\033[1;31m *** {msg} *** \033[0m\n')

    def parseString(self, hdfString):
        '''
        Reformat an hdf string
        '''
        return np.array(hdfString).item().decode()

    def parseVector(self, index, field='orbit'):
        '''
        Parse a vector from meta data with field and index

        Parameters
        ----------
        index : str
            keyword indcx.
        field : str optional
            field. The default is 'orbit'.

        Returns
        -------
        np array
            parsed vector/array.

        '''
        try:
            metadata = self.h5[self.product]['metadata']
            return np.array(metadata[field][index])
        except Exception:
            print(f'Could not parse science:{self.sar}:{self.product}:'
                  f'metadata:{field}:{index}')

    def parseDateStr(self, dateStr):
        '''
        Parse a date str

        Parameters
        ----------
        dateStr : TYPE
            DESCRIPTION.

        Returns
        -------
        dateAndTime : TYPE
            DESCRIPTION.
        date : TYPE
            DESCRIPTION.

        '''
        sLen = min(26, len(dateStr))
        dateStr = dateStr[0:sLen]
        # parse date
        try:
            dateAndTime = datetime.strptime(dateStr, '%Y-%m-%dT%H:%M:%S.%f')
            date = dateAndTime.strftime('%Y-%m-%d')
            return dateAndTime, date
        except Exception:
            self.printError('Could not parse date')
            return None, None

    def parseRefDate(self):
        '''
        Parse the reference date

        Returns
        -------
        None.

        '''
        image = 'reference'
        if self.isSecondary:
            image = 'secondary'
        #
        dateStr = self.parseString(
            self.h5['identification'][f'{image}ZeroDopplerStartTime'])
        #
        self.datetime, self.Date = self.parseDateStr(dateStr)

    def openHDF(self, hdfFile, referenceOrbitXML=None, secondaryOrbitXML=None):
        '''
        Open hdf and save self.h5Full and truncate to self.h5

        Parameters
        ----------
        hdfFile : str
            hdf file name.
        Returns
        -------
        None.

        '''
        if not os.path.exists(hdfFile):
            self.myerror(f'{hdfFile} does not exist')
        # Update XMLs
        if referenceOrbitXML is not None:
            self.referenceOrbitXML = referenceOrbitXML
        if secondaryOrbitXML is not None:
            self.secondaryOrbitXML = secondaryOrbitXML
        # Open hdf file
        self.hdfFile = hdfFile
        self.h5Full = h5py.File(hdfFile, 'r')
        # Truncate to remove tags common to all
        self.h5 = self.h5Full['science'][self.sar]
        # Set image name
        self.ImageName = os.path.basename(hdfFile)
        # Parse primpary parameters.
        self.parseParams()

    def getGranuleNames(self):
        '''
        Get the granule names

        Returns
        -------
        None.

        '''
        procInfo = self.h5[self.product]['metadata']['processingInformation']
        self.referenceGranule = \
            self.parseString(procInfo['inputs']['l1ReferenceSlcGranules'])
        self.secondaryGranule = \
            self.parseString(procInfo['inputs']['l1SecondarySlcGranules'])

    def getLookDirection(self):
        '''
        Get look direction

        Returns
        -------
        None.

        '''
        self.LookDirection = np.array(
            self.h5['identification']['lookDirection']).item().decode().lower()
        self.lookSign = {'right': 1.0, 'left': -1.0}[self.LookDirection]

    def getSlantRange(self):
        '''
        Input H5 and return a dicitionary with slant range spacing and
        near, middle, and far range.

        Returns
        -------
        None.

        '''
        # save bands to shorten lines
        bands = self.h5[self.product][self.bands]
        #
        if self.isSecondary:
            self.printError('warning secondary slant range not correct')
        # get the actual data
        self.slantRangeData = \
            bands[self.frequency][self.productType]['slantRange']
        # Near, middle and far
        self.MLNearRange = self.slantRangeData[0]
        self.MLCenterRange = (self.slantRangeData[0] +
                              self.slantRangeData[-1]) * 0.5
        self.MLFarRange = self.slantRangeData[-1]
        self.MLRangePixelSize = self.slantRangeData[1] - self.slantRangeData[0]

    def getZeroDopplerTimeSecondary(self):
        '''
        Get zero doppler time for the secondary image
        '''
        earlyTime, _ = self.parseDateStr(self.parseString(
            self.h5['identification']['secondaryZeroDopplerStartTime']))
        #
        # SLC first time
        zeroTime = {'hour': 0, 'minute': 0, 'second': 0, 'microsecond': 0}
        slcFirstZeroDopplerTime = (earlyTime - earlyTime.replace(**zeroTime)
                                   ).total_seconds()
        # convert to multi look time
        self.firstZeroDopplerTime = slcFirstZeroDopplerTime + \
            0.5 * (self.NumberAzimuthLooks - 1)/self.PRF
        #
        # Last zero doppler time
        lastTime, _ = self.parseDateStr(self.parseString(
            self.h5['identification']['secondaryZeroDopplerEndTime']))
        #
        # Use first earlyTime date
        slcLastZeroDopplerTime = (lastTime - earlyTime.replace(**zeroTime)
                                  ).total_seconds()
        #
        nSLCSamples = int((slcLastZeroDopplerTime -
                           slcFirstZeroDopplerTime) * self.PRF) + 1
        #
        self.MLAzimuthSize = int(nSLCSamples / self.NumberAzimuthLooks)
        #
        self.zeroDopplerTimeDelta = self.NumberAzimuthLooks / self.PRF
        #
        self.lastZeroDopplerTime = self.firstZeroDopplerTime + \
            (self.MLAzimuthSize - 1) * self.zeroDopplerTimeDelta
        self.midZeroDopplerTime = 0.5 * (self.firstZeroDopplerTime +
                                         self.lastZeroDopplerTime)
        #
        # compute the nominal time
        self.NominalTime = str(timedelta(
            seconds=np.around(self.firstZeroDopplerTime, decimals=5)))

    def getZeroDopplerTime(self, secondary=False):
        '''
        Input y5 and return dictionary with zero Doppler spacing and
        first, middle, and last zero Doppler.

        Returns
        -------
        None.

        '''
        bands = self.h5[self.product][self.bands]
        # Get the time array
        self.zeroDopplerTimeData = \
            bands[self.frequency][self.productType]['zeroDopplerTime']
        #
        # get start, mid, and end times
        self.firstZeroDopplerTime = self.zeroDopplerTimeData[0]
        self.midZeroDopplerTime = (self.zeroDopplerTimeData[0] +
                                   self.zeroDopplerTimeData[-1]) * 0.5
        self.lastZeroDopplerTime = self.zeroDopplerTimeData[-1]
        #
        # get the time interval per sample
        self.zeroDopplerTimeDelta = np.array(
            bands[self.frequency][self.productType]['zeroDopplerTimeSpacing']
            ).item()
        # compute the nominal time
        self.NominalTime = str(timedelta(
            seconds=np.around(self.firstZeroDopplerTime, decimals=5)))

    def getNumberOfLooks(self):
        '''
        Get number of looks in range and azimuth for multilook product.

        Returns
        -------
        None.

        '''
        # break up long dict
        metadata = self.h5[self.product]['metadata']
        parameters = metadata['processingInformation']['parameters']
        #
        self.NumberRangeLooks = np.array(
            parameters[self.productType][self.frequency]['numberOfRangeLooks']
            ).item()
        self.NumberAzimuthLooks = np.array(
            parameters[self.productType][self.frequency][
                'numberOfAzimuthLooks']
            ).item()

    def getSingleLookPixelSize(self):
        '''
        Get the single look pixel size

        Returns
        -------
        None.

        '''
        frequency = self.h5[self.product][self.bands][self.frequency]
        #
        self.MLAzimuthPixelSize = \
            np.array(frequency['sceneCenterAlongTrackSpacing']).item()
        #
        self.SLCAzimuthPixelSize = \
            self.MLAzimuthPixelSize / self.NumberAzimuthLooks
        # Range
        self.SLCRangePixelSize = self.MLRangePixelSize / self.NumberRangeLooks

    def getSize(self):
        '''
        Get number of looks in range and azimuth for multilook product.

        Returns
        -------
        None.

        '''
        frequency = self.h5[self.product][self.bands][self.frequency]
        polarization = frequency[self.productType][self.polarization]
        if self.layer is not None:
            polarization = polarization[self.layer]
        if not self.isSecondary:
            self.MLAzimuthSize, self.MLRangeSize = \
                polarization[self.productData].shape
        else:
            self.MLRangeSize = polarization[self.productData].shape[1]
            self.printError('Secondary size not complete')

    def getOrbitAndFrame(self):
        '''
        Get product frame and orbit.

        Returns
        -------
        None.

        '''
        #
        self.frame = np.array(self.h5['identification']['frameNumber']).item()
        self.orbit = \
            np.array(self.h5['identification']['absoluteOrbitNumber']).item()

    def getOffsetParams(self):
        '''
        Get offset params, a0, r0, deltaA, deltaR

        Returns
        -------
        None.

        '''
        procInfo = self.h5[self.product]['metadata']['processingInformation']
        frequency = procInfo['parameters'][self.productType][self.frequency]
        #
        setattr(self, 'deltaA',
                np.array(frequency['alongTrackSkipWindowSize']).item())
        setattr(self, 'deltaR',
                np.array(frequency['slantRangeSkipWindowSize']).item())
        setattr(self, 'a0',
                np.array(frequency['alongTrackStartPixel']).item())
        setattr(self, 'r0',
                np.array(frequency['slantRangeStartPixel']).item())

    def getOffsetLayerParams(self, layers=['layer1', 'layer2', 'layer3']):
        '''
        Get the window sizes for the individual layers

        Parameters
        ----------
        layers : TYPE, optional
            DESCRIPTION. The default is ['layer1', 'layer2', 'layer3'].

        Returns
        -------
        None.

        '''
        procInfo = self.h5[self.product]['metadata']['processingInformation']
        frequency = procInfo['parameters'][self.productType][self.frequency]
        #
        self.layerParams = {}
        #
        for layer in layers:
            self.layerParams[layer] = {}
            self.layerParams[layer]['alongTrackSearchWindowSize'] = \
                np.array(frequency[layer]['alongTrackSearchWindowSize']).item()
            self.layerParams[layer]['alongTrackWindowSize'] = \
                np.array(frequency[layer]['alongTrackWindowSize']).item()
            self.layerParams[layer]['slantRangeSearchWindowSize'] = \
                np.array(frequency[layer]['slantRangeSearchWindowSize']).item()
            self.layerParams[layer]['slantRangeWindowSize'] = \
                np.array(frequency[layer]['slantRangeWindowSize']).item()

    def effectivePRF(self, offsets=False):
        '''
        Compute an effective PRF for compatibility with other GrIMP software
            given a zero Doppler spacing and number of azimuth looks.

        Returns
        -------
        None.

        '''
        params = self.h5[self.product]['metadata']['processingInformation'][
            'parameters']
        if not self.isSecondary:
            singleLookDelta = np.array(
                params['reference']['frequencyA']['zeroDopplerTimeSpacing']
                ).item()
        else:
            singleLookDelta = np.array(
               params['secondary']['frequencyA']['zeroDopplerTimeSpacing']
               ).item()
        #
        # Compute the effective PRF
        self.PRF = 1. / singleLookDelta

    def lltoxy(self, lat, lon):
        '''
        Convert to lat/lon to xy
        '''
        if not hasattr(self, 'epsg'):
            self.getEPSG()
        if not hasattr(self, 'lltoxyXform'):
            self.lltoxyXform = pyproj.Transformer.from_crs("EPSG:4326",
                                                           f"EPSG:{self.epsg}")
        return self.lltoxyXform.transform(lat, lon)

    def _neighborDotProduct(self, x1, x2, x3, y1, y2, y3):
        '''
        Compute dot product of line segments centered on x2,y2 to find corner
        points.

        Parameters
        ----------
        x1, x2, x3: float
            x coords for points 1, 2, & 3.
        y1, y2, y3: float
            y coords for points 1, 2, & 3.

        Returns
        -------
        normDot : floatt
            Normalized do product ||(pt1-pt2) . (p3-pt2)||.

        '''
        dxa, dxb = x1-x2, x3-x2
        dya, dyb = y1-y2, y3-y2
        s1 = np.sqrt((dxa*dxa + dya*dya)*(dxb*dxb + dyb*dyb))
        normDot = np.abs((dxa*dxb + dya*dyb)/s1)
        return normDot

    def _getRectangle(self):
        '''
        From a bounding polygon from a NISAR HDF with vertices between corners,
        extract the polygon with only the corner vertices.

        Parameters
        ----------
        boundingPolygon : POLYGON
            Bounding polygon from hdf.

        Returns
        -------
        x, y: np float array
            corner points in xy coordinates.

        '''
        x, y = np.zeros(4), np.zeros(4)
        if self.epsg is None:
            self.getEPSG()

        # Compute corners
        i = 0
        for r in [self.MLNearRange, self.MLFarRange]:
            for t in [self.firstZeroDopplerTime, self.lastZeroDopplerTime]:
                x[i], y[i] = self.lltoxy(*self.RTtoLatLon(r, t, 0)[0:2])
                i += 1
        return x, y

    def _idCorners(self, x, y):
        '''
        Identify the corners as ll, lr, ur, ul and save in dict
        self.corners = {'ll': (lat, lon), 'lr'....}
        Parameters
        ----------
        x, y : float arrays
            x&y coordinate of corners in arbitrary order.

        Returns
        -------
        None.

        '''
        #
        # Get two most extreme points in each diretion.
        w = np.argsort(x)[0:2]
        e = np.flip(np.argsort(x))[0:2]
        s = np.argsort(y)[0:2]
        n = np.flip(np.argsort(y))[0:2]
        #
        # Find the points coorner points as intersection of two extreme points.
        ll = np.intersect1d(w, s)[0]
        lr = np.intersect1d(e, s)[0]
        ul = np.intersect1d(w, n)[0]
        ur = np.intersect1d(e, n)[0]
        #
        # Convert back to lat lon
        self.xytoll = pyproj.Transformer.from_crs(f"EPSG:{self.epsg}",
                                                  "EPSG:4326").transform
        lat, lon = self.xytoll(x, y)
        # save in dict
        self.corners = {'ll': (lat[ll], lon[ll]), 'lr': (lat[lr], lon[lr]),
                        'ur': (lat[ur], lon[ur]), 'ul': (lat[ul], lon[ul])}

    def getEPSG(self):
        '''
        Get epsg
        Returns
        -------
        Return epsg.
        '''
        self.epsg = np.array(
            self.h5[self.product]['metadata']['geolocationGrid']['epsg']
            ).item()

    def getCorners(self):
        '''
        Extract corners from bounding polygon with redudant points.
        '''
        # Get epsg
        self.getEPSG()
        #
        # Extract minimum rectangle.
        x, y = self._getRectangle()
        #
        # Sort and id corners
        self._idCorners(x, y)

    # def formatStateVectors(self):
    #     '''
    #     Write state vectors to dict with keys that match geodat format.

    #     Returns
    #     -------
    #     None.

    #     '''
    #     self.stateVectors = {}
    #     self.stateVectors['NumberOfStateVectors'] = self.NumberOfStateVectors
    #     self.stateVectors['TimeOfFirstStateVector'] = \
    #         self.TimeOfFirstStateVector
    #     self.stateVectors['StateVectorInterval'] = self.StateVectorInterval
    #     for i, pos, vel in zip(range(1, len(self.position)+1),
    #                            self.position, self.velocity):
    #         self.stateVectors[f'SV_Pos_{i}'] = pos
    #         self.stateVectors[f'SV_Vel_{i}'] = vel

    def getRangeErrorCorrection(self):
        '''
        get the range error correction, for now just 0

        Returns
        -------
        None.

        '''
        self.RangeErrorCorrection = 0.0

    def getTimeToFirstSample(self):
        '''
        For now this is zero and is just for compatability with geodat

        Returns
        -------
        None.

        '''
        self.TimeToFirstSLCSample = 0.0

    def getSkewOffsets(self):
        '''
        For now this is zero and is just for compatability with geodat

        Returns
        -------
        None.

        '''
        self.SkewOffset = 0.0

    def getSquint(self):
        '''
        For now this is zero and is just for compatability with geodat

        Returns
        -------
        None.

        '''
        self.Squint = 0.0

    def getDeltaT(self):
        '''
        Time correction sometimes used for other SARs, but for now is 0.

        Returns
        -------
        None.

        '''
        self.deltaT = 0.0

    def parseStateVectors(self, XMLOrbit=None):
        '''
        Setup an nisarOrbit instance to contain the state vectors

        Returns
        -------
        orbit  nisarOrbit.
            NISAR orbit instance with state vector information

        '''
        if XMLOrbit is None:
            if not self.isSecondary:
                h5OrbitGroup = self.h5[self.product]['metadata']['orbit']
            else:
                if self.debug:
                    self.printError('Using ref. orbit state vectors for debug')
                    h5OrbitGroup = self.h5[self.product]['metadata']['orbit']
                else:
                    self.printError('Add second orbit from h5 if implemented')
            print('oo', h5OrbitGroup)
            orbit = nisarOrbit(h5OrbitGroup=h5OrbitGroup)
        else:
            orbit = nisarOrbit(XMLOrbit=XMLOrbit)
        return orbit

    def getCorrectedTime(self):
        '''
        Correct time so the first time refers to the first slc pixel since
        GrIMP software uses SLC based times and corrects for ML internally

        Returns
        -------
        None.

        '''
        # Half ML pixel correction
        timeCorrection = 0.5 * (self.NumberAzimuthLooks - 1)/self.PRF
        correctedTime = self.firstZeroDopplerTime - timeCorrection
        correctedTimeString = str(timedelta(
            seconds=np.around(correctedTime, decimals=5)))
        # get pieces
        x = [float(x) for x in correctedTimeString.split(':')]
        # format as string
        self.CorrectedTime = f'{int(x[0]):02} {int(x[1]):02} {x[2]}'

    def getIncidenceAngleCube(self):
        '''
        Create an interpolator for the incidence angle data cube

        Returns
        -------
        None.

        '''
        metadata = self.h5[self.product]['metadata']
        incidenceAngleData = \
            np.array(metadata['geolocationGrid']['incidenceAngle'])
        slantRange = np.array(metadata['geolocationGrid']['slantRange'])
        zeroDopplerTime = \
            np.array(metadata['geolocationGrid']['zeroDopplerTime'])
        heightAboveEllipsoid = \
            np.array(metadata['geolocationGrid']['heightAboveEllipsoid'])
        #
        self.incidenceAngleCube = scipy.interpolate.RegularGridInterpolator(
                        (heightAboveEllipsoid, zeroDopplerTime, slantRange),
                        incidenceAngleData)

    def getElevationAngleCube(self):
        '''
        Create an interpolator for the incidence angle data cube

        Returns
        -------
        None.

        '''
        metadata = self.h5[self.product]['metadata']
        elevationAngleData = \
            np.array(metadata['geolocationGrid']['elevationAngle'])
        slantRange = np.array(metadata['geolocationGrid']['slantRange'])
        zeroDopplerTime = \
            np.array(metadata['geolocationGrid']['zeroDopplerTime'])
        heightAboveEllipsoid = \
            np.array(metadata['geolocationGrid']['heightAboveEllipsoid'])
        #
        self.elevationAngleCube = scipy.interpolate.RegularGridInterpolator(
                        (heightAboveEllipsoid, zeroDopplerTime, slantRange),
                        elevationAngleData)

    def computeAngles(self, slantRange, zeroDopplerTime, zWGS84,
                      degrees=False):
        '''
        Compute elevation angle for slantRange, zeroDoppler, zWGS84.

        Parameters
        ----------
        slantRange : float
            slant range.
        zeroDopplerTime : floate
            zero Doppler time.
        z : float
            Elevation (WGS84 ellipsoid) of point .
        degrees : TYPE, optional
            Return values degrees. The default is False.

        Returns
        -------
        elevationAngle : float
            Elevation angle for point.
        incidenceAngle : float
            Incidence angle for point.

        '''

        lat, lon, z = self.RTtoLatLon(slantRange,  zeroDopplerTime, zWGS84)
        position, velocity = self.getSatPositionAndVel(zeroDopplerTime)
        ReH = np.linalg.norm(position)
        #
        if not hasattr(self, 'LLTOECEF'):
            self.LLtoECEF = pyproj.Transformer.from_crs("EPSG:4326",
                                                        "EPSG:4978",).transform
        #
        ptRadius = np.linalg.norm(self.LLtoECEF(lat, lon, z))
        #
        elevationAngle = np.arccos((slantRange**2 + ReH**2 - ptRadius**2) /
                                   (2.0 * slantRange * ReH))
        incidenceAngle = np.arcsin(np.sin(elevationAngle) * ReH/ptRadius)
        #
        if degrees:
            return np.degrees(elevationAngle), np.degrees(incidenceAngle)
        return elevationAngle, incidenceAngle

    def getCenterIncidenceAngle(self):
        '''
        Get center incidence angle

        Returns
        -------
        None.

        '''
        self.MLElevationAngle, self.MLIncidenceCenter = self.computeAngles(
            self.MLCenterRange, self.midZeroDopplerTime, 0, degrees=True)

    def getXYCube(self):
        '''
        Create interpolators for x and y coordinates x or y = f(h, zdt, sr)

        Returns
        -------
        None.

        '''
        # Get the coordinate data
        metadata = self.h5[self.product]['metadata']
        coordinateXData = np.array(metadata['geolocationGrid']['coordinateX'])
        coordinateYData = np.array(metadata['geolocationGrid']['coordinateY'])
        # Get the slant range, time, and height indices
        slantRange = np.array(metadata['geolocationGrid']['slantRange'])
        zeroDopplerTime = \
            np.array(metadata['geolocationGrid']['zeroDopplerTime'])
        heightAboveEllipsoid = \
            np.array(metadata['geolocationGrid']['heightAboveEllipsoid'])
        #
        # Compute the the interpolators
        self.coordinateXCube = scipy.interpolate.RegularGridInterpolator(
                        (heightAboveEllipsoid, zeroDopplerTime, slantRange),
                        coordinateXData)
        self.coordinateYCube = scipy.interpolate.RegularGridInterpolator(
                        (heightAboveEllipsoid, zeroDopplerTime, slantRange),
                        coordinateYData)

    def getCenterLatLon(self):
        '''
        Get the center lat/lon for image (at elevation =0)

        Returns
        -------
        None.

        '''
        self.CenterLatLon = self.RTtoLatLon(self.MLCenterRange,
                                            self.midZeroDopplerTime, 0)[0:2]

    def getSatelliteHeight(self):
        '''
        Compute the nominal spacecraft altitude mid scene

        Returns
        -------
        None.

        '''
        if not hasattr(self, 'ECEFtoLL'):
            self.ECEFtoLL = pyproj.Transformer.from_crs("EPSG:4978",
                                                        "EPSG:4326").transform
        # Get the ECEF coords
        x = self.orbit.xsv([self.midZeroDopplerTime])[0]
        y = self.orbit.ysv([self.midZeroDopplerTime])[0]
        z = self.orbit.zsv([self.midZeroDopplerTime])[0]
        # Convert to geodetic lat, lon, height above ellipsoid
        latSat, lonSat, hSat = self.ECEFtoLL(x, y, z)
        self.SpaceCraftAltitude = hSat

    def getOrbitPassDirection(self):
        '''
        Extract the orbit pass direction (Ascending/Descending)

        Returns
        -------
        None.

        '''
        passType = self.parseString(
            self.h5['identification']['orbitPassDirection'])
        self.PassType = {'ASCEND': 'Ascending',
                         'DESCEND': 'Descending'
                         }[passType]

    def getWavelength(self):
        '''
        Get the wavelength

        Returns
        -------
        None.

        '''
        bands = self.h5[self.product][self.bands]
        # get frequency
        self.centerFrequency = \
            np.array(bands[self.frequency]['centerFrequency']).item()
        # And compute the wavelength
        self.Wavelength = self.cLight / self.centerFrequency

    def genGeodatProperties(self):
        '''
        Create dictionary with all of the properties in a dict ready
        for geojson output

        Returns
        -------
        None.

        '''
        self.geodatDict = {}
        # These are tags for the geojson properties
        keys = ['ImageName', 'Date', 'NominalTime', 'NumberRangeLooks',
                'NumberAzimuthLooks',
                'MLRangeSize', 'MLAzimuthSize', 'PRF', 'MLNearRange',
                'MLCenterRange', 'MLFarRange', 'RangeErrorCorrection',
                'LookDirection', 'PassType', 'CenterLatLon',
                'TimeToFirstSLCSample', 'SkewOffset', 'Squint',
                'EarthRadiusMajor', 'EarthRadiusMinor', 'MLIncidenceCenter',
                'SpaceCraftAltitude', 'CorrectedTime', 'Wavelength',
                'SLCRangePixelSize', 'SLCAzimuthPixelSize', 'deltaT']
        #
        # Get values for these keys
        for key in keys:
            self.geodatDict[key] = getattr(self, key)
        #
        # Now append the state vectors
        for sv in self.orbit.stateVectors:
            svData = self.orbit.stateVectors[sv]
            # Convert to list
            if type(svData) is np.ndarray:
                svData = list(svData)
            self.geodatDict[sv] = svData

    def formatGeojson(self, myString):
        '''
        Format a geojson string for readability.

        Parameters
        ----------
        myString : str
            GeoJson unformatted string

        Returns
        -------
        formatted geoJson string

        '''
        # indent function
        def myIndent(n, x, formattedString):
            for space in [' '] * n * 2:
                formattedString += space
            return formattedString
        #
        # Start bracket counds at 0
        braceCount, bracketCount = 0, 0
        # Remote space
        myString = myString.replace(' \"', '\"')
        #
        # Start with empty string
        formattedString = ''
        # And add characters back with additional formatting as neeed.
        for x in myString:
            # Add return and indent if not in list (leftBracket==0)
            if x in ',' and bracketCount == 0:
                formattedString += x
                formattedString += '\n'
                formattedString = myIndent(braceCount, x, formattedString)
            # Update bracket count if left bracket
            elif x == '[':
                bracketCount += 1
                formattedString += x
            # decrement left bracket count
            elif x == ']':
                bracketCount -= 1
                formattedString += x
            # Update brace count if left brace and add return and indent
            elif x == '{':
                braceCount += 1
                formattedString += x
                formattedString += '\n'
                formattedString = myIndent(braceCount, x, formattedString)
            # Decrement braceCount and add new line and indent
            elif x == '}':
                formattedString += '\n'
                braceCount -= 1
                formattedString = myIndent(braceCount, x, formattedString)
                formattedString += x
            # Regulary character so append
            else:
                formattedString += x
        return formattedString

    def writeGeodatGeojson(self, filename=None, path='.', secondary=False):
        '''
        Write a geodat geojson file

        Parameters
        ----------
        filename : str, optional
            filename. The default is None, which will create
            geodatNRxNA.geojson.
        path : str, optional
            Path where file will be created. The default is '.'.

        Returns
        -------
        None.

        '''
        if not hasattr(self, 'geodatDict'):
            self.genGeodatProperties()
        if not hasattr(self, 'corners'):
            self.getCorners()
        #
        geoJsonGeometry = geojson.Polygon(
            [[self.corners[x] for x in ['ll', 'lr', 'ur', 'ul', 'll']]])
        self.geodatGeojson = geojson.Feature(geometry=geoJsonGeometry,
                                             properties=self.geodatDict)
        #
        if filename is None:
            filename = f'geodat{self.NumberRangeLooks}x' \
                    f'{self.NumberAzimuthLooks}.geojson'
        if secondary and hasattr(self, 'secondary'):
            filenameSecondary = filename.replace('.geojson',
                                                 '.secondary.geojson')
            self.secondary.writeGeodatGeojson(filename=filenameSecondary)
        #
        # Write and format the string, then write to the file
        geojsonString = self.formatGeojson(geojson.dumps(self.geodatGeojson))
        with open(os.path.join(path, filename), 'w') as fpGeojson:
            print(geojsonString, file=fpGeojson)

    def getImageData(self, productField, layer=None):
        '''

        Parameters
        ----------
        productField : str
            Product field to extract (eg. unwrappedPhase).
        layer : str, optional
            Layer key if needed (mostly for offsets). The default is None.

        Returns
        -------
        None.

        '''
        # get bands group
        bands = self.h5[self.product][self.bands]
        data = \
            bands[self.frequency][self.productType][self.polarization]
        # Apply layer if needed.
        if layer is not None:
            data = data[layer]
        # set data value
        setattr(self, productField, np.array(data[productField]))

    def _writeVrt(self, newVRTFile, sourceFiles, descriptions,
                  byteOrder=None, bands=0, eType=gdal.GDT_Float32,
                  geoTransform=[0.5, 0.5, 1., 0., 0., 1.], metaData=None,
                  setSRS=False, noDataValue=-2.0e9):
        '''
        Write a vrt for the file. Note sourcefiles and descriptions have
        to be passed in.
        '''
        if type(sourceFiles) is not list:
            sourceFiles = [sourceFiles]
        if type(descriptions) is not list:
            descriptions = [descriptions]
        # get the meta data from the template
        # Kill any old file
        if os.path.exists(newVRTFile):
            os.remove(newVRTFile)
        # Create VRT
        drv = gdal.GetDriverByName("VRT")
        vrt = drv.Create(newVRTFile, self.MLRangeSize, self.MLAzimuthSize,
                         bands=bands, eType=eType)

        vrt.SetGeoTransform(geoTransform)
        #
        if setSRS:
            sr = osr.SpatialReference()
            sr.ImportFromEPSG(self.epsg)
            vrt.SetSpatialRef(sr)

        if metaData is None:
            metaData = {}
        if byteOrder is None:
            if "ByteOrder" in metaData:
                byteOrder = metaData["ByteOrder"]
            else:
                byteOrder = "MSB"
                metaData["ByteOrder"] = byteOrder
        else:
            metaData["ByteOrder"] = byteOrder
        if metaData is not None:
            vrt.SetMetadata(metaData)
        # Look to add bands
        for sourceFile, description, bandNumber in \
                zip(sourceFiles, descriptions, range(1, 1 + len(sourceFiles))):
            options = [f"SourceFilename={sourceFile}", "relativeToVRT=1",
                       "subclass=VRTRawRasterBand", f"BYTEORDER={byteOrder}",
                       bytes(0)]
            vrt.AddBand(eType, options=options)
            band = vrt.GetRasterBand(bandNumber)
            band.SetMetadataItem("Description", description)
            band.SetNoDataValue(noDataValue)
        # Close the vrt
        vrt = None

    def _writeImageData(self, fileName, x, dataType, noData=-2.e9):
        '''
        write a binary image of size nx by ny with dataType = to one
        ['f4','>f4','>u2','u2','>i2','i2','>u4','u4','>i4','i4','u1']
        '''
        #
        # binary types
        types = ['f4', '>f4', 'f8', '>f8', '>u2', 'u2', '>i2', 'i2', '>u4',
                 'u4', '>i4', 'i4', 'u1']
        if dataType not in types:
            self.printError(f'writeImage: invalid data type - {dataType}')
        #
        # Type translation
        types = {'f4': 'float32', 'f8': 'float64', 'i2': 'int16',
                 'i4': 'int32', 'u1': 'uint8', 'u2': 'uint16', 'u4': 'uint32'}
        x1 = x.astype(types[dataType.replace('>', '')])
        x1[np.isnan(x1)] = noData
        if '>' in dataType:
            x1.byteswap(True)

        with open(fileName, 'w') as fpOut:
            x1.tofile(fpOut)

    def getGeoCoordinates(self):
        '''
        Get Coordinate information (origin, size, and spacing) and asssemble
        information for a geodat file.

        Returns
        -------
        None.

        '''
        productType = \
            self.h5[self.product][self.bands][self.frequency][self.productType]
        # Epsg
        self.epsg = np.array(
            self.h5[self.product][self.bands][self.frequency]['projection']
            ).item()
        # coordinate setup
        self.xCoords = np.array(productType['xCoordinates'])
        self.yCoords = np.array(productType['yCoordinates'])
        self.x0 = np.min(self.xCoords)
        self.y0 = np.min(self.yCoords)
        self.nx = len(self.xCoords)
        self.ny = len(self.yCoords)
        self.dx = np.array(productType['xCoordinateSpacing']).item()
        self.dy = np.array(productType['yCoordinateSpacing']).item()
        self.geodat = [[self.nx, self.ny],
                       [np.abs(self.dx), np.abs(self.dy)],
                       [self.x0/1000, self.y0/1000.]]
        print(self.geodat)
        print(self.epsg)
        print(self.x0, self.y0, self.nx, self.ny, self.dx, self.dy)

    def writeGeodat(self, filename):
        '''
        Write geolocation information to a GrIMP geodat file.
        '''
        with open(filename, 'w') as fpOut:
            print('# 2', file=fpOut)
            print(';\n; Image size (pixels) nx ny\n;', file=fpOut)
            print(f'{self.nx:d} {self.ny:d}', file=fpOut)
            print(';\n; Pixel size (m) deltaX deltaY\n;', file=fpOut)
            print(f'{np.abs(self.dx):.4f} {np.abs(self.dy):.4f}', file=fpOut)
            print(';\n; Origin, lower left corner (km) Xo  Yo\n;', file=fpOut)
            print('{self.x0/1000.:.4f} {self.y0/1000.:.4f}', file=fpOut)
            print('&', file=fpOut)

    def _computeLookError(self, elook, ph, pv, position, slantRange,
                          ReMajorZ2, ReMinorZ2):
        '''
        Compute look angle error for lat/lon determination

        Parameters
        ----------
        elook : float
            current look angle.
        ph : 3 element vector
            horizontal component of the unit pointing vector.
        pv : 3 element vector
            vertical component of the unit pointing vector..
        position : 3 element vector
            satellite position
        slantRange : float
            Slant range.
        ReMajorZ2 : float
            (ReMajor + z)**2.
        ReMinorZ2 : float
            (ReMajor + z)**2.
        Returns
        -------
        delook : float
            Look angle error.
        t : 3 element vector
            target location.
        dt : 3 element vector
            Differential position.

        '''
        #
        # Latest pointing vector
        p = np.sin(elook) * ph + np.cos(elook) * pv
        #
        # Target location using latest look pointing vector
        t = position + slantRange * p
        txy, tz = t[0:2], t[2]
        #
        # Error to minimize
        error = 1.0 - np.dot(txy, txy) / ReMajorZ2 - (tz**2) / ReMinorZ2
        return error, t

    def SMLocateZD(self, position, velocity, slantRange, z):
        '''
        Based on old Soren Madsen algorithm to compute lat/lon from slant
        range and altitude given sat positoin

        Parameters
        ----------
        position : 3 element vector
            Satellite ECEF coordinates.
        velocity : 3 element vector
            Satellite velocity.
        slantRange : float
            Slant range.
        z : float
            WGS84 elevation

        Returns
        -------
        lat, lon, z
            geocentric coordinates

        '''
        #
        # Calculate the vertical unit vector
        Rsat = np.linalg.norm(position)
        nhat = position / Rsat
        #
        # Calculate the horizontal component of the unit pointing vector
        # which is orthogonal to the velocity vector
        ph = np.cross(velocity, nhat)
        phhat = ph / np.linalg.norm(ph)
        #
        # Calculate the vertical component of the unit pointing vector
        # which is perpendicular to the velocity vector
        pv = np.cross(velocity, phhat)
        pvhat = pv / np.linalg.norm(pv)
        #
        #  Calculate effective earth radius and satellite altitude
        posxy, posz = position[0:2], position[2]
        reff = np.sqrt(Rsat**2 /
                       (np.dot(posxy, posxy) / self.EarthRadiusMajor**2 +
                        (posz**2) / self.EarthRadiusMinor**2))
        #
        # Compute the approximate look angle using reff at target loc
        temp = (Rsat**2.0 + slantRange**2 - (reff + z)**2) / \
            (2.0 * Rsat * slantRange)
        elook = self.lookSign * np.arccos(temp)
        #
        # Minimize look error function
        ReMajorZ2 = (self.EarthRadiusMajor + z)**2
        ReMinorZ2 = (self.EarthRadiusMinor + z)**2
        def lookOptimize(elook, phhat, pvhat, position, slantRange,
                         ReMajorZ2, ReMinorZ2):
            return self._computeLookError(elook, phhat, pvhat, position,
                                          slantRange, ReMajorZ2, ReMinorZ2)[0]
        #
        elook = optimize.leastsq(lookOptimize, elook,
                                 args=(phhat, pvhat, position, slantRange,
                                       ReMajorZ2, ReMinorZ2))[0][0]
        #
        # Target location using latest final pointing vector
        # compute final values
        p = np.sin(elook) * phhat + np.cos(elook) * pvhat
        t = position + slantRange * p
        #
        return self.ECEFtoLL(*list(t))

    def getSatPositionAndVel(self, zdTime):
        '''
        Get satellite position for given zdTime

        Parameters
        ----------
        zdTime : float
            Zero Doppler time

        Returns
        -------
        None.

        '''
        position = [self.orbit.xsv([zdTime])[0],
                    self.orbit.ysv([zdTime])[0],
                    self.orbit.zsv([zdTime])[0]]
        velocity = [self.orbit.vxsv([zdTime])[0],
                    self.orbit.vysv([zdTime])[0],
                    self.orbit.vzsv([zdTime])[0]]
        return position, velocity

    def RTtoLatLon(self, slantRange, zdTime, z):
        '''
        Basic R/D coords to lat lon at a point.

        Parameters
        ----------
        slantRange : float
            Slant range in m.
        zdTime : float
            zero Doppler time in s.
        z : float
            WGS84 elevation.

        Returns
        -------
        lat, lon, z
            Lat/lon/z coordinates.

        '''
        if not hasattr(self, 'ECEFtoLL'):
            self.ECEFtoLL = pyproj.Transformer.from_crs("EPSG:4978",
                                                        "EPSG:4326").transform
        #
        # User position/velocity state vector interpolation for time
        position, velocity = self.getSatPositionAndVel(zdTime)
        # Return lat/lon/z for slantRange, z
        return self.SMLocateZD(position, velocity, slantRange, z)
