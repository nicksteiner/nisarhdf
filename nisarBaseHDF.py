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
import shapely.wkt
import shapely
import geojson
import pyproj
from osgeo import gdal, osr
import os


class nisarBaseHDF():
    '''
    Abstract class to define parser for NISAR HDF products.
    '''

    __metaclass__ = ABCMeta

    def __init__(self, sar='LSAR', product='RUNW', frequency='frequencyA',
                 productType='interferogram', polarization='HH', layer=None,
                 productData='unwrappedPhase', bands='swaths'):
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
        self.initTags(sar=sar,
                      product=product,
                      frequency=frequency,
                      productType=productType,
                      polarization=polarization,
                      layer=layer,
                      productData=productData,
                      bands=bands)
        # Constants
        self.EarthRadiusMajor = 6378137.0
        self.EarthRadiusMinor = 6356752.3142
        self.f = (self.EarthRadiusMajor -
                  self.EarthRadiusMinor) / self.EarthRadiusMajor
        self.cLight = 2.99792458e8

    @abstractmethod
    def parseParams():
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

    def geodeticLL(self, x, y, z):
        '''
        Parameters
        ----------
        x, y, z : float
            ECF coordinates.

        Returns
        -------
        lat, lon, h: float
            geodatic lat/lon/h.
        '''
        p = np.sqrt(x**2 + y**2)
        q = np.arctan2(z * self.EarthRadiusMajor, p * self.EarthRadiusMinor)
        esq = self.f * (2 - self.f)
        eps = esq / (1 - esq)
        print(self.f, esq, eps)
        sinq3, cosq3 = np.sin(q)**3, np.cos(q)**3

        phi = np.arctan2(z + eps * self.EarthRadiusMinor * sinq3,
                         p - esq * self.EarthRadiusMajor * cosq3)
        # Comput lat/lon
        lat = np.rad2deg(phi)
        lon = np.rad2deg(np.arctan2(y, x))
        # Compute height
        v = self.EarthRadiusMajor / np.sqrt(1 - esq * np.sin(phi)**2)
        h = (p / np.cos(phi)) - v
        return lat, lon, h

    def parseRefDate(self):
        '''
        Parse the reference date

        Returns
        -------
        None.

        '''
        dateStr = self.parseString(
            self.h5['identification']['referenceZeroDopplerStartTime'])
        # parse date
        try:
            self.datetime = datetime.strptime(dateStr, '%Y-%m-%dT%H:%M:%S.%f')
            self.Date = self.datetime.strftime('%Y-%m-%d')
        except Exception:
            self.printError('Could not parse date')

    def openHDF(self, hdfFile):
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
        self.hdfFile = hdfFile
        self.h5Full = h5py.File(hdfFile, 'r')
        # Truncate to remove tags common to all
        self.h5 = self.h5Full['science'][self.sar]
        self.parseParams()

    def getLookDirection(self):
        '''
        Get look direction

        Returns
        -------
        None.

        '''
        self.LookDirection = np.array(
            self.h5['identification']['lookDirection']).item().decode()

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
        # get the actual data
        self.slantRangeData = \
            bands[self.frequency][self.productType]['slantRange']
        # Near, middle and far
        self.MLNearRange = self.slantRangeData[0]
        self.MLCenterRange = (self.slantRangeData[0] +
                              self.slantRangeData[-1]) * 0.5
        self.MLFarRange = self.slantRangeData[-1]
        self.MLRangePixelSize = self.slantRangeData[1] - self.slantRangeData[0]

    def getZeroDopplerTime(self):
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
        self.MLAzimuthSize, self.MLRangeSize = \
            polarization[self.productData].shape

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
        # hack, remove when fields are populated.
        #if self.frame == 0:

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
        # Make sure necessary params defined
        if not hasattr(self, 'zeroDopplerTimeDelta'):
            self.getZeroDopplerTime()
        if offsets:
            if not hasattr(self, 'deltaA'):
                self.getOffsetParams()
            singleLookDelta = self.zeroDopplerTimeDelta / self.deltaA
        else:
            if not hasattr(self, 'NumberAzimuthLooks'):
                self.getNumberOfLooks()
            singleLookDelta = \
                self.zeroDopplerTimeDelta / self.NumberAzimuthLooks
        #
        # Compute the effective PRF

        self.PRF = 1. / singleLookDelta

    def llToImageXY(self, lat, lon):
        '''
        Convert to lat/lon to xy
        '''
        lltoxy = pyproj.Transformer.from_crs("EPSG:4326", f"EPSG:{self.epsg}")
        x, y = lltoxy.transform(lat, lon)
        return x, y

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

    def _getRectangle(self, boundingPolygon):
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
        x, y = [], []
        tol = .1
        tolR = 500
        print(f'\n\t\033[1;31m Applying tolerance of {tol} s and {tolR} m to '
              'keep in geolocation grid, remove for fixed products \033[0m\n')
        if not hasattr(self, 'coordinateXCube'):
            self.getXYCube()
        for r in [self.MLNearRange+tolR, self.MLFarRange-tolR]:
            for t in [self.firstZeroDopplerTime + tol,
                      self.lastZeroDopplerTime - tol]:
                x.append(self.coordinateXCube([0, t, r])[0])
                y.append(self.coordinateYCube([0, t, r])[0])
        #
        # last point repeated twice in sample product
        # lonP, latP = np.array(boundingPolygon.exterior.xy)[:, 0:-1]
        # xx, yy = self.llToImageXY(latP, lonP)
        # print(xx, yy)
        # #
        # for i in range(0, len(xx)-1):
        #     normDot = self._neighborDotProduct(xx[i-1], xx[i], xx[i+1],
        #                                        yy[i-1], yy[i], yy[i+1])
        #     if normDot < .5:
        #         x.append(xx[i])
        #         y.append(yy[i])
        return np.array(x), np.array(y)

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

    def getCorners(self):
        '''
        Extract corners from bounding polygon with redudant points.
        '''
        # Extract polygon
        boundingPolygon = shapely.wkt.loads(np.array(
            self.h5['identification']['boundingPolygon']).item().decode())
        # Get epsg
        metadata = self.h5[self.product]['metadata']
        self.epsg = np.array(metadata['geolocationGrid']['epsg']).item()
        #
        # Extract minimum rectangle.
        x, y = self._getRectangle(boundingPolygon)
        #
        # Sort and id corners
        self._idCorners(x, y)

    def formatStateVectors(self):
        '''
        Write state vectors to dict with keys that match geodat format.

        Returns
        -------
        None.

        '''
        self.stateVectors = {}
        self.stateVectors['NumberOfStateVectors'] = self.NumberOfStateVectors
        self.stateVectors['TimeOfFirstStateVector'] = \
            self.TimeOfFirstStateVector
        self.stateVectors['StateVectorInterval'] = self.StateVectorInterval
        for i, pos, vel in zip(range(1, len(self.position)+1),
                               self.position, self.velocity):
            self.stateVectors[f'SV_Pos_{i}'] = pos
            self.stateVectors[f'SV_Vel_{i}'] = vel

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

    def parseStateVectors(self):
        '''
        Parse the state vectors. Store with names compatiable with geodat
        files.

        Returns
        -------
        None.

        '''
        stateVectors = {}
        # Check that state vectors don't start the day before
        if not hasattr(self, 'datetime'):
            self.parseRefDate(self)
        stateVectors['date'] = self.datetime
        # Get time array
        time = self.parseVector('time', field='orbit')
        #
        # use it to compute other time related values
        self.NumberOfStateVectors = len(time)
        self.TimeOfFirstStateVector = time[0]
        self.StateVectorInterval = time[1] - time[0]
        #
        # Error check: compute time from params and confirm matches time array
        t = [self.TimeOfFirstStateVector +
             self.StateVectorInterval * x for x in range(0, len(time))]
        if np.std(t-time) > 0.0001:
            self.printError('Non uniform state vector intervals')
        #
        # Get the orbit type
        metadata = self.h5[self.product]['metadata']
        self.stateVectorType = self.parseString(metadata['orbit']['orbitType'])
        #
        # Get position and velocity
        self.position = self.parseVector('position', field='orbit')
        self.velocity = self.parseVector('velocity', field='orbit')
        #
        # Create interpolationrs pos and vel for state vectors
        for i, pos, vel in zip(range(0, 3),
                               ['xsv', 'ysv', 'zsv'],
                               ['vxsv', 'vysv', 'vzsv']):
            setattr(self, pos,
                    scipy.interpolate.RegularGridInterpolator(
                        [time], self.position[:, i], method='quintic'))
            setattr(self, vel,
                    scipy.interpolate.RegularGridInterpolator(
                        [time], self.velocity[:, i], method='quintic'))
        # done
        self.formatStateVectors()

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

    def getCenterIncidenceAngle(self):
        '''
        Get center incidence angle

        Returns
        -------
        None.

        '''
        if not hasattr(self, 'incidenceAngleCube'):
            self.getIncidenceAngleCube()
        self.MLIncidenceCenter = self.incidenceAngleCube(
            [0, self.midZeroDopplerTime, self.MLCenterRange])[0]

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
        if not hasattr(self, 'coordinateXCube'):
            self.getXYCube()
        #
        # Get image center in x, y coords
        xc = self.coordinateXCube(
            [0, self.midZeroDopplerTime, self.MLCenterRange])[0]
        yc = self.coordinateYCube(
            [0, self.midZeroDopplerTime, self.MLCenterRange])[0]
        #
        # Make sure xytoll xform exists
        if not hasattr(self, 'xytoll'):
            self.getCorners()
        #
        self.CenterLatLon = self.xytoll(xc, yc)

    def getSatelliteHeight(self):
        '''
        Compute the nominal spacecraft altitude mid scene

        Returns
        -------
        None.

        '''
        # Get the ECEF coords
        x = self.xsv([self.midZeroDopplerTime])[0]
        y = self.ysv([self.midZeroDopplerTime])[0]
        z = self.zsv([self.midZeroDopplerTime])[0]
        # Convert to geodetic lat, lon, height above ellipsoid
        latSat, lonSat, hSat = self.geodeticLL(x, y, z)
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
        keys = ['Date', 'NominalTime', 'NumberRangeLooks',
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
        for sv in self.stateVectors:
            svData = self.stateVectors[sv]
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

    def writeGeodatGeojson(self, filename=None, path='.'):
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
