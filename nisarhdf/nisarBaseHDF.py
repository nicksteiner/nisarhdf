#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 08:47:03 2024

@author: ian
"""
from abc import ABCMeta, abstractmethod
import h5py
import numpy as np
import pyproj
from osgeo import gdal, osr
from scipy import optimize
import os
from nisarhdf.nisarOrbit import nisarOrbit
import geopandas as gpd
import rasterio
from datetime import datetime
from nisarhdf.writeMultiBandVrt import writeMultiBandVrt
import gc

gdal.UseExceptions()

gdalTypes = {'float32': gdal.GDT_Float32,
             'float64': gdal.GDT_Float64,
             'complex64':  gdal.GDT_CFloat32,
             'uint16': gdal.GDT_UInt16,
             'int8': gdal.GDT_Byte,
             'uint8': gdal.GDT_Byte}


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
        self.initTags(**keywords)
        # Constants WGS84
        self.EarthRadiusMajor = 6378137.0
        self.EarthRadiusMinor = 6356752.3142
        self.EarthF = 1./298.257223563
        self.f = (self.EarthRadiusMajor -
                  self.EarthRadiusMinor) / self.EarthRadiusMajor
        self.cLight = 2.99792458e8
        self.epsg = None
        #
        self.commonParams = ['referenceOrbit', 'datetime', 'referenceGranule',
                             'secondaryOrbit', 'secondaryDatetime',
                             'secondaryGranule', 'frame', 'LookDirection',
                             'PassType', 'Wavelength', 'PRF', 'epsg']
        #
        self.noDataValuesTiff = {'connectedComponents': 65535,
                                 'coherenceMagnitude': np.nan,
                                 'ionospherePhaseScreen': np.nan,
                                 'ionospherePhaseScreenUncertainty': np.nan,
                                 'unwrappedPhase': np.nan,
                                 'wrappedInterferogram': np.nan + 1j*np.nan,
                                 'slantRangeOffset': np.nan,
                                 'slantRangeOffsetVariance': np.nan,
                                 'alongTrackOffset': np.nan,
                                 'alongTrackOffsetVariance': np.nan,
                                 'crossOffsetVariance': np.nan,
                                 'correlationSurfacePeak': np.nan,
                                 'snr': np.nan}

        self.noDataValuesBinary = {'connectedComponents': 65535,
                                   'coherenceMagnitude': -2e9,
                                   'ionospherePhaseScreen': -2e9,
                                   'ionospherePhaseScreenUncertainty': -2e9,
                                   'unwrappedPhase': -2e9,
                                   'wrappedInterferogram': -2e9 + 1j * -2e9,
                                   'slantRangeOffset': -2e9,
                                   'slantRangeOffsetVariance': -2e9,
                                   'alongTrackOffset': -2e9,
                                   'alongTrackOffsetVariance': -2e9,
                                   'crossOffsetVariance': -2e9,
                                   'correlationSurfacePeak': -2e9,
                                   'snr': -2e9}

    #
    # Abtract method placeholders
    #
    @abstractmethod
    def parseParams(**keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass

    @abstractmethod
    def setupDataCube(self, *args, **keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass

    @abstractmethod
    def interpGrid(self, *args, **keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass

    @abstractmethod
    def getEPSG(**keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass

    @abstractmethod
    def getExtent(**keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass

    @abstractmethod
    def getGeoTransform(**keywords):
        ''' Abstract pointPlot method - define in child class. '''
        pass
    #
    # Various "utility" methods for printing, type conversion, coordinate
    # conversion etc.
    #

    def refresh(self, dataFields=None):
        '''
        Reload data to undo any modifications

        Parameters
        ----------
        dataFields : list, optional
            List of fields. The default is None, which forces all.

        Returns
        -------
        None.

        '''
        if dataFields is None:
            dataFields = self.dataFields
        if self.product in ['ROFF', 'GOFF']:
            self.getLayers(dataFields)
        else:
            self.loadData(dataFields)

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

    def ensureNP(self, x):
        '''
        Make sure variable is an np array

        Parameters
        ----------
        x : list or nd.array
            Convert list to numpy array if needed.

        Returns
        -------
        None.

        '''
        if isinstance(x, list):
            x = np.array(x)
        return(x)

    def assembleMeta(self, bands=None):
        '''
        Create a dictionary with all of the meta data to save

        Returns
        -------
        None.

        '''
        self.meta = {}
        for param in self.commonParams + self.levelParams + self.productParams:
            self.meta[param] = getattr(self, param)
        if bands is not None:
            self.meta['bands'] = bands

    def printParams(self):
        '''
        Print some of the more important parameters

        Returns
        -------
        None.

        '''
        # Print parameters
        for param in self.commonParams + self.levelParams + self.productParams:
            print(f'{self.product}.{param} = {getattr(self, param)} ')
        # List data fields
        print('\nData Fields')
        for field in self.dataFields:
            print(f'{self.product}.{field}')

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

    def toScalar(self, hdfValue):
        '''
        Convert rscalar hdf value to regular value
        '''
        return np.array(hdfValue).item()

    def parseDateStr(self, dateStr):
        '''
        Parse a date str

        Parameters
        ----------
        dateStr : str
            Date str.

        Returns
        -------
        dateAndTime : datetime
            date and time as datetime.
        date : str
            Date (without time) as string.

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
        #
        if not self.isSecondary:
            dateStr = self.parseString(
                self.h5['identification']['secondaryZeroDopplerStartTime'])
            self.secondaryDatetime, self.secondaryDate = self.parseDateStr(
                dateStr)

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

    def xytoll(self, lat, lon):
        '''
        Convert to lat/lon to xy
        '''
        if not hasattr(self, 'epsg'):
            self.getEPSG()
        if not hasattr(self, 'xytollXform'):
            self.xytollXform = pyproj.Transformer.from_crs(f"EPSG:{self.epsg}",
                                                           "EPSG:4326")
        return self.xytollXform.transform(lat, lon)

    #
    # Open hdf and extract parameters.
    #

    def openHDF(self, hdfFile, referenceOrbitXML=None, secondaryOrbitXML=None,
                referenceOrbit=None, secondaryOrbit=None, **keywords):
        '''
        Open hdf and save self.h5Full and truncate to self.h5

        Parameters
        ----------
        hdfFile : str
            hdf file.
        referenceOrbitXML : str, optional
            XML orbit file to override HDF statevectors. The default is None.
        secondaryOrbitXML : str, optional
            XML orbit file to override HDF statevectors.. The default is None.
        referenceOrbit : int, optional
            Override secondary orbit in HDF (used in development where missing
            in hdf). The default is None.
        secondaryOrbit : int, optional
            Override secondary orbit in HDF (used in development where missing
            in hdf). The default is None.
        **keywords : TYPE
            keywords to pass to parse params.

        Returns
        -------
        None.

        '''
        if not os.path.exists(hdfFile):
            self.myerror(f'{hdfFile} does not exist')
        # Update XMLs
        for attr, value in zip(['referenceOrbitXML', 'secondaryOrbitXML'],
                               [referenceOrbitXML, secondaryOrbitXML]):
            if value is not None:
                setattr(self, attr, value)
            else:
                setattr(self, attr, None)
        # Open hdf file
        self.hdfFile = hdfFile
        self.h5Full = h5py.File(hdfFile, 'r')
        # Truncate to remove tags common to all
        self.h5 = self.h5Full['science'][self.sar]
        # Set image name
        self.ImageName = os.path.basename(hdfFile)
        # Parse primary parameters.
        self.parseParams(referenceOrbit=referenceOrbit,
                         secondaryOrbit=secondaryOrbit,
                         **keywords)

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
        self.LookDirection = self.parseString(
            self.h5['identification']['lookDirection']).lower()
        self.lookSign = {'right': 1.0, 'left': -1.0}[self.LookDirection]

    def getNumberOfLooks(self):
        '''
        Get number of looks in range and azimuth for multilook product.

        Returns
        -------
        None.

        '''
        # break up long dict
        metadata = self.h5[self.product]['metadata']
        productType =\
            metadata['processingInformation']['parameters'][self.productType]
        #
        self.NumberRangeLooks = self.toScalar(
            productType[self.frequency]['numberOfRangeLooks'])
        self.NumberAzimuthLooks = self.toScalar(
             productType[self.frequency]['numberOfAzimuthLooks'])

    def getOrbitAndFrame(self, referenceOrbit=None,
                         secondaryOrbit=None, frame=None):
        '''
        Get product frame and orbit.

        Returns
        -------
        None.

        '''
        #
        if frame is None:
            self.frame = self.toScalar(
                self.h5['identification']['frameNumber'])
        else:
            self.frame = frame
        #
        if referenceOrbit is None:
            self.referenceOrbit = \
                self.toScalar(
                    self.h5['identification']['absoluteOrbitNumber'])
        else:
            self.referenceOrbit = referenceOrbit
        #
        if secondaryOrbit is None:
            print("Cannot sread secondary orbit yet")
        else:
            self.secondaryOrbit = secondaryOrbit

    def effectivePRF(self):
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
            singleLookDelta = self.toScalar(
                params['reference']['frequencyA']['zeroDopplerTimeSpacing'])
        else:
            singleLookDelta = self.toScalar(
               params['secondary']['frequencyA']['zeroDopplerTimeSpacing'])
        #
        # Compute the effective PRF
        self.PRF = 1. / singleLookDelta

    def getSLCSlantRange(self):
        '''
        Get SLC slantRange information

        Returns
        -------
        None.

        '''
        # save bands to shorten lines
        metadata = self.h5[self.product]['metadata']
        parameters = metadata['processingInformation']['parameters']
        #
        if not self.isSecondary:
            imageParams = parameters['reference']['frequencyA']
        else:
            imageParams = parameters['secondary']['frequencyA']
        self.SLCRangeSize = self.toScalar(imageParams['numberOfRangeSamples'])
        #
        # Get near range and SLC spacing for either reference or seconeary
        self.SLCNearRange = self.toScalar(imageParams['slantRangeStart'])
        self.SLCRangePixelSize = self.toScalar(
            imageParams['slantRangeSpacing'])
        self.SLCNumberOfRangeSamples = self.toScalar(
            imageParams['numberOfRangeSamples'])
        self.SLCFarRange = self.SLCNearRange + \
            (self.SLCRangeSize - 1) * self.SLCRangePixelSize
        self.SLCCenterRange = (self.SLCNearRange + self.SLCFarRange) * 0.5

    def getSLCZeroDopplerTime(self, secondary=False):
        '''
        Get Zero doppler start time from the SLC information in the metdaa

        Parameters
        ----------
        isSecondary : booulan, optional
            indicates Secondary. The default is False.

        Returns
        -------
        None.

        '''
        zeroTime = {'hour': 0, 'minute': 0, 'second': 0, 'microsecond': 0}
        metaData = self.h5[self.product]['metadata']
        parameters = metaData['processingInformation']['parameters']
        if secondary:
            SLC = 'secondary'
            self.SLCAzimuthSize = self.toScalar(
                parameters['secondary']['frequencyA']['numberOfAzimuthLines'])
        else:
            SLC = 'reference'
            self.SLCAzimuthSize = self.toScalar(
                parameters['reference']['frequencyA']['numberOfAzimuthLines'])
        #
        earlyTime, _ = self.parseDateStr(self.parseString(
            parameters[SLC][self.frequency]['zeroDopplerStartTime']))
        #
        # SLC first time
        self.SLCFirstZeroDopplerTime = \
            (earlyTime - earlyTime.replace(**zeroTime)).total_seconds()
        self.SLCZeroDopplerTimeSpacing = self.toScalar(
            parameters[SLC][self.frequency]['zeroDopplerTimeSpacing'])
        self.SLCLastZeroDopplerTime = self.SLCFirstZeroDopplerTime + \
            (self.SLCAzimuthSize - 1) * self.SLCZeroDopplerTimeSpacing
        self.SLCMidZeroDopplerTime = (self.SLCFirstZeroDopplerTime +
                                      self.SLCLastZeroDopplerTime) * 0.5
        # Compute scene center along track spacing
        self.SLCSceneCenterAlongTrackSpacing()

    def SLCSceneCenterAlongTrackSpacing(self):
        '''
        Compute self.SLCAzimuthPixelSize directly since sample data products
        had errors and this value was not provided for all data types.

        Returns
        -------
        None.

        '''
        if not hasattr(self, 'orbit'):
            self.orbit = self.parseStateVectors()
        # Compute positions as scent center delta ZDT apart
        p1 = self.RTtoECEF(self.SLCCenterRange,
                           self.SLCMidZeroDopplerTime,
                           0)
        p2 = self.RTtoECEF(self.SLCCenterRange,
                           self.SLCMidZeroDopplerTime +
                           self.SLCZeroDopplerTimeSpacing,
                           0)
        self.SLCAzimuthPixelSize = np.sqrt(np.sum((p1-p2)**2))

    def getOffsetWindowParams(self, layers=['layer1', 'layer2', 'layer3']):
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
        self.windowParams = {}
        #
        for layer in layers:
            self.windowParams[layer] = {}
            self.windowParams[layer]['alongTrackSearchWindowSize'] = \
                np.array(frequency[layer]['alongTrackSearchWindowSize']).item()
            self.windowParams[layer]['alongTrackWindowSize'] = \
                np.array(frequency[layer]['alongTrackWindowSize']).item()
            self.windowParams[layer]['slantRangeSearchWindowSize'] = \
                np.array(frequency[layer]['slantRangeSearchWindowSize']).item()
            self.windowParams[layer]['slantRangeWindowSize'] = \
                np.array(frequency[layer]['slantRangeWindowSize']).item()

    def getOrbitPassDirection(self):
        '''
        Extract the orbit pass direction (Ascending/Descending)

        Returns
        -------
        None.

        '''
        passType = self.parseString(
            self.h5['identification']['orbitPassDirection']).lower()
        # print(passType)
        self.PassType = {'ascending': 'ascending',
                         'descending': 'descending'
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
        self.centerFrequency = self.toScalar(
            bands[self.frequency]['centerFrequency'])
        # And compute the wavelength
        self.Wavelength = self.cLight / self.centerFrequency

    def parseStateVectors(self, XMLOrbit=None, SLC=False):
        '''
        Setup an nisarOrbit instance to contain the state vectors

        Returns
        -------
        orbit  nisarOrbit.
            NISAR orbit instance with state vector information

        '''
        if XMLOrbit is None:
            if not self.isSecondary:
                h5OrbitGroup = \
                    self.h5[self.product]['metadata']['orbit']
                if not SLC:
                    h5OrbitGroup = h5OrbitGroup['reference']
            else:
                h5OrbitGroup = \
                    self.h5[self.product]['metadata']['orbit']['secondary']
            orbit = nisarOrbit(
                h5OrbitGroup=h5OrbitGroup,
                firstZeroDopplerTime=self.SLCFirstZeroDopplerTime,
                lastZeroDopplerTime=self.SLCLastZeroDopplerTime)
        else:
            orbit = nisarOrbit(
                XMLOrbit=XMLOrbit,
                firstZeroDopplerTime=self.SLCFirstZeroDopplerTime,
                lastZeroDopplerTime=self.SLCLastZeroDopplerTime)
            print('orbit parsed')
        return orbit

    #
    # IO stuff, extract data from hdf, read, and write files.
    #

    def _writeMatchFile(self, filenameRoot, data, layer, matchType):
        '''
        Write file with match type for each pixel (GrIMP only)

        Parameters
        ----------
        fileNameRoot : TYPE
            DESCRIPTION.
        data : TYPE
            DESCRIPTION.
        matchType : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        mt = np.zeros(np.squeeze(data[layer-1, :, :]).shape, 'byte')
        mt[np.squeeze(data[layer-1, :, :]) > -2.e8] = matchType
        self._writeImageData(f'{filenameRoot}.layer{layer}.mt', mt, 'u1')
        # write the mt vrt
        descriptions = [f'matchType {layer}']
        writeMultiBandVrt(f'{filenameRoot}.layer{layer}.mt.vrt',
                          self.OffsetRangeSize, self.OffsetAzimuthSize,
                          [f'{filenameRoot}.layer{layer}.mt'],
                          descriptions,
                          eType=gdal.GDT_Byte,
                          geoTransform=[-0.5, 1., 0., -0.5, 0., 1.],
                          noDataValue=0)

    def getProducts(self, productFields):
        '''
        Load the product specificed product fields for offsets if not already
        done.

        Parameters
        ----------
        productFields : list of str
            Product fields to load.

        Returns
        -------
        None.
       '''

        for productField in productFields:
            # Get layers if need
            if not hasattr(self, productField):
                print(f'getting layers for {productField}')
                self.getLayers([productField])

    def getImageData(self, productField, layer=None):
        '''
        Get image data corresponding to productField and create copy as
        self.productField.

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

    def getLayers(self, layerTypes, layers=['layer1', 'layer2', 'layer3']):
        '''
        Get Layers for a given set layerTypes (e.g, ['alongTrackOffset' ...]).

        Save as self.layerType[numberOfLayers, :, :]

        Parameters
        ----------
        layerTypes : list of str
            Names of the layer variables to extract.
        layers : list, optional
            The layers to extract. The default is ['layer1', 'layer2',
                                                   'layer3'].

        Returns
        -------
        None.

        '''
        h5Data = \
            self.h5[self.product][self.bands][self.frequency][self.productType]
        self.dataFields = []
        for layerType in layerTypes:
            layerData = []
            self.dataFields.append(layerType)
            for layer in layers:
                layerData.append(h5Data[self.polarization][layer][layerType])
            setattr(self, layerType, np.stack(layerData))
        tmp = np.squeeze(getattr(self, layerTypes[0])[len(layers) - 1, :, :])
        self.noDataLocations = np.isnan(tmp)

    def loadData(self, fields):
        '''
        Load data to np arrays.

        Returns
        -------
        None.

        '''
        self.dataFields = []
        for field in fields:
            self.getImageData(field)
            self.dataFields.append(field)
        self.noDataLocations = np.isnan(getattr(self, fields[0]))

    def findNoDataValue(self, band, tiff):
        '''
        Determine no data value from either predefined values or default
        to np.nan (tiff=True) or -2e9 (tiff=False)

        Parameters
        ----------
        band : stra
            band name.
        tiff : bool
            True for tiff file.

        Returns
        -------
        number of appropriate type
            no data value.
        '''
        if tiff:
            if band in self.noDataValuesTiff:
                return self.noDataValuesTiff[band]
            else:
                return np.nan  # Default if not defined
        else:
            if band in self.noDataValuesBinary:
                return self.noDataValuesBinary[band]
            else:
                return -2.e9

    def readGeodatAsDataFrame(self, filename=None):
        '''
        Read Geodat files a pandas data frme

        Parameters
        ----------
        filename : str, optional
            Name of geodat file. The default is None, which then uses
            geodatNumberRangeLooksxNumberAzimuthLooks.geojson

        Returns
        -------
        None.

        '''
        if filename is None:
            filename = f'geodat{self.NumberRangeLooks}x' \
                    f'{self.NumberAzimuthLooks}.geojson'
        #
        self.geoDF = gpd.read_file(filename)

    def writeData(self, filenameRoot, bands=None, tiff=True,
                  byteOrder='LSB', grimp=False,
                  suffixes=None, layers=[1, 2, 3],
                  matchTypes=[1, 2, 3], saveMatch=False,
                  scaleToPixels=False,  geojsonName=None,
                  geojsonNameSecondary=None):
        '''
        Write data to binary or tiff file for all data types. Non offset
        results are saved as individual files (filenameRoot.band[.tif] which
        can be read via a single vrt file (filenameRoot.vrt). If save as tiffs,
        .tif is appended.

        For offsets, files are saved as filenameRoot.layer.suffix[.tif], where
        suffixes are determined by the band to suffix dictionary:
        {'slantRangeOffset': '.dr','slantRangeOffsetVariance': '.vr',
         'alongTrackOffset': '.da', 'alongTrackOffsetVariance': '.va',
         'crossOffsetVariance': '.vc', 'correlationSurfacePeak': '.cc',
         'snr': '.snr'}.

        For either offsets or non-offset data, interpolated fields are saved
        as filenameRoot.fieldName[.tif] (e.g. filenameRoot.incidenceAngle)

        Parameters
        ----------
        filenameRoot : str
            Each band will be save as filenameRoot.band(.tif) along with
            filenameRoot.vrt.
        bands : list of str, optional
            Bands to save (e.g., ['unwrappedPhase]. The default is None,
            which will dump all of the data fields.
        tiff : bool, optional
            Write data to tiff files. The default is True.
        byteOrder : str, optional
            Allow binary file to be written with 'LSB' or 'MSB'. The default
            is 'LSB'.
        grimp : bool, optional
            For GrIMP compatability, force binary files with geocoded data to
            use a lower-left corner geotransform. The default is False.

        Offset Only
        ------------
        suffixes : list of str, optional
            suffixes for output offset fields . The default is ['.da', '.va',
            '.dr', '.va', '.cc', '.snr'].
        layers : list of ints, optional
            Layers to save. The default is [1, 2, 3].
        saveMatch : boolean, optional
            Save a binary match type (mt) file for GrIMP workflow. The
            default is False.
        scaleToPixels : boolean, optional
            Scale offsets to pixel values for GrIMP workflow. The
            default is False.
        geojsonName : str optional
            geojson for redference image for GrIMP workflow. The default is
            None.
        geojsonNameSecondary : str, optional
            geojson for secondary image for GrIMP workflow. The default is
            None.
        Returns
        -------
        None.

        '''
        # This section makes sure that interpolated fields are written as
        # regular data while offsets with layers are treated specially.
        isData, isDerived = False, False
        for band in bands:
            if band in self.dataFields:
                isData = True
            else:
                isDerived = True
        if isDerived and isData and self.product in ['ROFF', 'GOFF']:
            self.printError('Cannot write offset data and interplated results'
                            ' at the same time')
            return
        # Now write with appropriate routine
        if self.product not in ['ROFF', 'GOFF'] or isDerived:
            self._writeNonOffsetData(filenameRoot,
                                     bands=bands,
                                     tiff=tiff,
                                     grimp=grimp)
        else:
            self._writeOffsets(filenameRoot,
                               bands=bands,
                               tiff=tiff,
                               byteOrder=byteOrder,
                               grimp=grimp,
                               suffixes=suffixes,
                               layers=layers,
                               matchTypes=matchTypes,
                               saveMatch=saveMatch,
                               scaleToPixels=scaleToPixels,
                               geojsonName=geojsonName,
                               geojsonNameSecondary=geojsonNameSecondary)

    def _writeOffsets(self, filenameRoot,
                      bands=None,
                      tiff=True,
                      byteOrder='LSB',
                      suffixes=None,
                      grimp=False,
                      layers=[1, 2, 3],
                      matchTypes=[1, 2, 3],
                      saveMatch=False,
                      scaleToPixels=False,
                      geojsonName=None,
                      geojsonNameSecondary=None):
        '''
        Save offsets as individual layer files with vrt for each layer that
        combines all of the bands.

        Parameters
        ----------
        filenameRoot : str
            Basename for offsets (e.g., offsets.dr, .da...).

        bands : list of str, optional
            Product Fields to save. The default is all.
        tiff : bool, optional
            Write data to tiff files. The default is True.
        byteOrder : str, optional
            Allow binary file to be written with 'LSB' or 'MSB'. The default
            is 'LSB'.
        grimp : bool, optional
            For GrIMP compatability, force binary files with geocoded data to
            use a lower-left corner geotransform. The default is False.
        suffixes : list of str, optional
            suffixes for output offset fields . The default is ['.da', '.va',
            '.dr', '.va', '.cc', '.snr'].
        layers : list of ints, optional
            Layers to save. The default is [1, 2, 3].
        saveMatch : boolean, optional
            Save a binary match type (mt) file for GrIMP workflow. The
            default is False.
        scaleToPixels : boolean, optional
            Scale offsets to pixel values for GrIMP workflow. The
            default is False.
        geojsonName : str optional
            geojson for redference image for GrIMP workflow. The default is
            None.
        geojsonNameSecondary : str, optional
            geojson for secondary image for GrIMP workflow. The default is
            None.
        Returns
        -------
        None.

        '''
        # Setup meta data
        if bands is None:
            bands = self.dataFields
        self.assembleMeta()
        meta = self.meta.copy()
        meta['bands'] = bands
        meta['ByteOrder'] = byteOrder
        if grimp:
            meta['geo1'] = geojsonName
            meta['geo2'] = geojsonNameSecondary
        # Make sure tiff appended to vrt files.
        tiffSuffix = {True: '.tif', False: ''}[tiff]
        #
        # Loop over layers with and create a vrt for each layer
        for layer, matchType in zip(layers, matchTypes):
            sourceFiles, descriptions, dataTypes, noDataValues = [], [], [], []
            #
            # Lo
            for band in bands:
                # If one of bands, use predefined suffies otherwise use band
                if band in self.suffixes:
                    suffix = self.suffixes[band]
                else:
                    suffix = band
                # Filename
                sourceFile = f'{filenameRoot}.layer{layer}{suffix}{tiffSuffix}'
                sourceFiles.append(sourceFile)
                descriptions.append(f'{band}')
                # condition data
                data = getattr(self, band)
                # Scale from meters to pixels if requested
                scaleFactor = 1.
                if scaleToPixels and band in self.scaleFactors:
                    scaleFactor = self.scaleFactors[band]
                bandData = np.squeeze(data[layer-1, :, :]) * scaleFactor
                dataTypes.append(gdalTypes[str(bandData.dtype)])
                noDataValues.append(self.findNoDataValue(band, tiff))
                self._writeImageData(sourceFile, bandData,
                                     tiff=tiff,
                                     noDataValue=noDataValues[-1],
                                     dataType=dataTypes[-1],
                                     byteOrder=byteOrder,
                                     grimp=grimp)
                #
                # If along track offsets save match type
                if band == 'alongTrackOffset' and saveMatch:
                    self._writeMatchFile(filenameRoot, bandData, layer,
                                         matchType)
            # Save the vrt
            sy, sx = bandData.shape
            epsg = None
            if self.coordType == 'GEO':
                epsg = self.epsg
            writeMultiBandVrt(f'{filenameRoot}.layer{layer}.vrt',
                              sx, sy,
                              sourceFiles,
                              descriptions,
                              eType=dataTypes,
                              geoTransform=self.getGeoTransform(grimp=grimp,
                                                                tiff=tiff),
                              noDataValue=noDataValues, metaData=meta,
                              byteOrder=byteOrder, tiff=tiff, epsg=epsg)
            #
            if grimp and not tiff:
                self.writeOffsetsDatFile(f'{filenameRoot}.layer{layer}.dat',
                                         geodat1=geojsonName,
                                         geodat2=geojsonNameSecondary)

    def _writeNonOffsetData(self, filenameRoot, bands=None, tiff=True,
                            byteOrder='LSB', grimp=False):
        '''
        Write non-offset data to binary or tiff file.

        Parameters
        ----------
        filenameRoot : str
            Each band will be save as filenameRoot.band(.tif) along with
            filenameRoot.vrt.
        bands : list of str, optional
            Bands to save (e.g., ['unwrappedPhase]. The default is None,
            which will dump all of the data fields.
        tiff : bool, optional
            Write data to tiff files. The default is True.
        byteOrder : str, optional
            Allow binary file to be written with 'LSB' or 'MSB'. The default
            is 'LSB'.
        grimp : bool, optional
            For GrIMP compatability, force binary files with geocoded data to
            use a lower-left corner geotransform. The default is False.

        Returns
        -------
        None.

        '''
        #
        if tiff and byteOrder != 'LSB':
            self.printError(f'Byte order {byteOrder} not supported for tiffs')
            return
        # Add '.tif' suffix to diff files
        suffix = {True: '.tif', False: ''}[tiff]
        # Default to all fields as bands
        if bands is None:
            bands = self.dataFields
        # Refresh meta data
        self.assembleMeta()
        meta = self.meta.copy()
        meta['bands'] = bands
        meta['ByteOrder'] = byteOrder
        # Results used to create vrt
        sourceFiles, descriptions, dataTypes, noDataValues = [], [], [], []
        # Loop and write each field
        for band in bands:
            filename = f'{filenameRoot}.{band}{suffix}'
            data = getattr(self, band)
            # Save bands nam, file, and data type
            descriptions.append(band)
            sourceFiles.append(filename)
            dataTypes.append(gdalTypes[str(data.dtype)])
            # Get no data value, generally np.nan for tiff and -2e9 for binary
            noDataValue = self.findNoDataValue(band, tiff)
            noDataValues.append(noDataValue)
            # Write the data to a binary or tiff file.
            self._writeImageData(filename, data,
                                 tiff=tiff,
                                 noDataValue=noDataValue,
                                 dataType=str(data.dtype),
                                 byteOrder=byteOrder,
                                 grimp=grimp)
        # Write a vrt
        sy, sx = data.shape
        writeMultiBandVrt(f'{filenameRoot}.vrt',
                          sx, sy,
                          sourceFiles,
                          descriptions,
                          eType=dataTypes,
                          geoTransform=self.getGeoTransform(grimp=grimp,
                                                            tiff=tiff),
                          noDataValue=noDataValues, metaData=meta,
                          byteOrder=byteOrder, tiff=tiff, epsg=self.epsg)

    def _writeImageData(self, filename, data,  byteOrder="LSB",
                        dataType=None, grimp=False, tiff=True,
                        noDataValue=None):
        '''
        Call routine to either write data as a tiff or flat binary file.

        Parameters
        ----------
        filename : TYPE
            File name for the data.
        data : np.array()
            Data to save to filename.
        byteOrder : str, optional
            Allow binary file to be written with 'LSB' or 'MSB'. The default
            is 'LSB'.
        grimp : bool, optional
            For GrIMP compatability, force binary files with geocoded data to
            use a lower-left corner geotransform. The default is False.
        tiff : bool, optional
            Write to tiff if True, otherwise flat binary. The default is True.
        noDataValue : numpber, optional
            The no data vaalue. The default is None.

        Returns
        -------
        None.

        '''
        if not tiff:
            self._writeBinaryImageData(filename, data,
                                       byteOrder=byteOrder,
                                       dataType=None,
                                       noDataValue=noDataValue,
                                       grimp=grimp)
            # This write geolocation info file for GrIMP workflows
            if grimp and self.coordType == 'GEO':
                self.writeGeodat(f'{filename}.geodat')
        else:
            self._writeTiffImageData(filename, data,
                                     dataType=None,
                                     noDataValue=noDataValue)

    def _writeTiffImageData(self, filename, data, dataType=None,
                            noDataValue=-2.e9):
        '''
        Write data as a tiff  file.

        Parameters
        ----------
        filename : TYPE
            File name for the data.
        data : np.array()
            Data to save to filename.
        dataType : str, optional
            Data type (e.g., 'float32'). The default is none, which will
            determine the type directly from the data.
        noDataValue : number, optional
            The no data vaalue. The default is None.

        Returns
        -------
        None.

        '''
        if dataType is None:
            dataType = str(data.dtype)
        # Tiff options
        options = ['BIGTIFF=NO', 'COMPRESS=LZW', 'GEOTIFF_VERSION=1.1',
                   'RESAMPLING=AVERAGE']
        if dataType != 'complex64':
            options.append('PREDICTOR=YES')
        # To create COG, first create a version in memory.
        driver = gdal.GetDriverByName('MEM')
        nRows, nColumns = data.shape
        dst_ds = driver.Create('', nColumns, nRows, 1, gdalTypes[dataType])
        dst_ds.SetGeoTransform(self.getGeoTransform(tiff=True))
        # set projection
        if self.coordType == 'GEO':
            sr = osr.SpatialReference()
            sr.ImportFromEPSG(self.epsg)
            dst_ds.SetProjection(sr.ExportToWkt())
        else:
            data = np.flipud(data)
        #
        # driver specific stuff to initalize band
        dst_ds.FlushCache()
        dst_ds.GetRasterBand(1).WriteArray(data)
        dst_ds.GetRasterBand(1).SetNoDataValue(noDataValue)
        #
        # Create copy for the COG.
        driver = gdal.GetDriverByName('COG')
        dst_ds2 = driver.CreateCopy(filename, dst_ds, options=options)
        dst_ds2.FlushCache()
        dst_ds.FlushCache()
        del(dst_ds)
        del(dst_ds2)
        gc.collect()

    def _writeBinaryImageData(self, filename, data, dataType=None,
                              byteOrder="LSB", noDataValue=-2.e9, grimp=False):
        '''
        write data as a flat binary image with specfied data type.

        Parameters
        ----------
        filename : str
            File name of file data is written to.
        data : np.array()
            The single band of data to save.
        dataType : str, optional
            Data type (e.g., 'float32'). The default is none, which will
            determine the type directly from the data.
        byteOrder : str, optional
            Allow binary file to be written with 'LSB' or 'MSB'. The default
            is 'LSB'.
        noDataValue : TYPE, optional
            DESCRIPTION. The default is -2.e9.
        grimp : bool, optional
            For GrIMP compatability, force binary files with geocoded data to
            use a lower-left corner geotransform so flipud before writing.
            The default is False.
        Returns
        -------
        None.

        '''
        if dataType is None:
            dataType = str(data.dtype)
        #
        x1 = data.astype(dataType)
        # replace nans with noDataValue
        x1[np.isnan(x1)] = noDataValue
        # Swap byte order if MSB
        if byteOrder == "MSB":
            print('swapping')
            x1.byteswap(True)
        # Flip geocoded images for grimp
        if grimp and self.coordType == 'GEO':
            x1 = np.flipud(x1)
            self.writeGeodat(f'{filename}.geodat')
        # Save data
        with open(filename, 'w') as fpOut:
            x1.tofile(fpOut)
    #
    # Discard outliers, apply masks and corrections.
    #

    def removeOutlierOffsets(self, filterField,
                             thresholds=[0.07, 0.05, 0.025],
                             layers=[1, 2, 3]):
        '''
        Remove outliners based on correlation or snr threshold

        Parameters
        ----------
        filterField : str
            Field to use to discard outliners ('correlationSurfacePeak' or
                                               'snr'.
        thresholds : TYPE, optional
            Thresholds use to discard outlier.
            The default is [0.07, 0.05, 0.025],
            which are defaults for correlation.
        layers : list, optional
            List of layers to process. The default is [1, 2, 3].

        Returns
        -------
        None.

        '''
        if self.product not in ['ROFF', 'GOFF']:
            print('fProduct type is {self.product} but Method only works with'
                  'ROFF and GOFF')
            return
        # Load the product fields if not already done
        # Loop over layers
        qaData = getattr(self, filterField)
        azOffsets = getattr(self, 'alongTrackOffset')
        rgOffsets = getattr(self, 'slantRangeOffset')
        for layer in layers:
            qaLayer = qaData[layer-1]
            bad = qaLayer < thresholds[layer-1]
            azOffsets[layer-1][bad] = np.nan
            rgOffsets[layer-1][bad] = np.nan

    def maskPhase(self, largest=True):
        '''
        Mask phase using connected components mask so all zero components
        are discarded.
        largest : Bool
            Retain only the largest connected component.
        Returns
        -------
        None.

        '''
        if self.product not in ['RUNW', 'GUNW']:
            print(f'Product type is {self.product} but Method only works with'
                  'ROFF and GOFF')
            return
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
            labels = labels[np.logical_and(labels > 0, labels < 65535)]
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

    def ionosphereCorrectPhase(self, masked=False):
        '''
        Correct phase for ionosphere

        Parameters
        ----------
        masked : bool, optional
            Apply the correction to the masked phase. The default is False.

        Returns
        -------
        None.

        '''
        if masked:
            self.ionosphereCorrectedPhase = \
                self.maskedUnwrappedPhase - self.ionospherePhaseScreen
        else:
            self.ionosphereCorrectedPhase = \
                self.unwrappedPhase - self.ionospherePhaseScreen
#
# Date cubes
#

    def incidenceAngleCube(self, x, y, z, maskNoData=None, save=False):
        '''
        Interpolate incidenceAngleCube data cube
        Parameters
        ----------
        x or slantTrange : nparray
            x (G) or slant (R) coordinates.
        y : nparray
            y (G) or slant (ZD) coordinates.
        z : nparray
            Elevation coordinates.
        save: bool, optional
            If true, save the result as self.result (e.g.,
            self.incidenceAngle). The default is False.
        Returns
        -------
        incidenceAngleCube : nparrays
            incidenceAngleCube at desired points

        '''
        if not hasattr(self, 'incidenceAngleInterpolator'):
            self.setupDataCube('incidenceAngle', zGridPts=None)
        incidenceAngle = self.interpGrid(self.incidenceAngleInterpolator,
                                         x, y, z, maskNoData=maskNoData)
        if save:
            self.incidenceAngle = incidenceAngle.astype(np.float32)
            self.incidenceAngle[np.isnan(self.incidenceAngle)] = np.nan

        return incidenceAngle

    def elevationAngleCube(self, x, y, z, maskNoData=None, save=False):
        '''
        Interpolate elevationAngle data cube
        Parameters
        ----------
        x or slantTrange : nparray
            x (G) or slant (R) coordinates.
        y : nparray
            y (G) or slant (ZD) coordinates.
        z : nparray
            Elevation coordinates.
        save: bool, optional
            If true, save the result as self.result (e.g.,
            self.elevationAngle). The default is False.
        Returns
        -------
        elevationAngle : nparrays
            elevationAngle at desired points

        '''
        if not hasattr(self, 'elevationAngleInterpolator'):
            self.setupDataCube('elevationAngle', zGridPts=None)
        elevationAngle = self.interpGrid(self.elevationAngleInterpolator,
                                         x, y, z, maskNoData=maskNoData)
        if save:
            self.elevationAngle = elevationAngle.astype(np.float32)
        return elevationAngle

    def groundTrackVelocityCube(self, x, y, z, maskNoData=None, save=False):
        '''
        Interpolate groundTrackVelocity data cube
        Parameters
        ----------
        x or slantTrange : nparray
            x (G) or slant (R) coordinates.
        y : nparray
            y (G) or slant (ZD) coordinates.
        z : nparray
            Elevation coordinates.
        save: bool, optional
            If true, save the result as self.result (e.g.,
            self.groundTrackVelocity). The default is False.
        Returns
        -------
        groundTrackVelocity : nparrays
            groundTrackVelocity at desired points

        '''
        if not hasattr(self, 'groundTrackVelocityInterpolator'):
            self.setupDataCube('groundTrackVelocity', zGridPts=None)
        groundTrackVelocity = self.interpGrid(
            self.groundTrackVelocityInterpolator, x, y, z,
            maskNoData=maskNoData)
        if save:
            self.groundTrackVelocity = groundTrackVelocity
        return groundTrackVelocity.astype(np.float32)

    def baselineCube(self, x, y, z, maskNoData=None, save=False):
        '''
        Interpolate baseline data cube

        Parameters
        ----------
        x or slantTrange : nparray
            x (G) or slant (R) coordinates.
        y : nparray
            y (G) or slant (ZD) coordinates.
        z : nparray
            Elevation coordinates.
        save: bool, optional
            If true, save the result as self.result (e.g., self.bn and
            self.bp). The default is False.
        Returns
        -------
        bn, bp : nparrays
            perpendicular and parallel components of baseline

        '''
        # Setup interpolators if needed.
        if not hasattr(self, 'perpendicularBaselineInterpolator'):
            self.setupDataCube('perpendicularBaseline', zGridPts=[0, -1])
        if not hasattr(self, 'parallelBaselineInterpolator'):
            self.setupDataCube('parallelBaseline', zGridPts=[0, -1])
        # interpolate
        bn = self.interpGrid(self.perpendicularBaselineInterpolator, x, y, z,
                             maskNoData=maskNoData)
        bp = self.interpGrid(self.parallelBaselineInterpolator, x, y, z,
                             maskNoData=maskNoData)
        if save:
            self.bn = bn.astype(np.float32)
            self.bp = bp.astype(np.float32)
        return bn, bp

    def alongTrackUnitVectorCube(self, x, y, z, maskNoData=None, save=False):
        '''
        Interpolate alongTrackUnitVector data cube

        Parameters
        ----------
        x or slantTrange : nparray
            x (G) or slant (R) coordinates.
        y : nparray
            y (G) or slant (ZD) coordinates.
        z : nparray
            Elevatio coordinates.
        save: bool, optional
            If true, save the result as self.result (e.g.,
            self.alongTrackUnitVectorX and self.alongTrackUnitVectorY). The
            default is False.
        Returns
        -------
        ATx, ATY : nparrays
            x and y components of the along track unit vector

        '''
        # Setup interpolators if needed.
        if not hasattr(self, 'alongTrackUnitVectorXInterpolator'):
            self.setupDataCube('alongTrackUnitVectorX', zGridPts=None)
        if not hasattr(self, 'alongTrackUnitVectorYInterpolator'):
            self.setupDataCube('alongTrackUnitVectorY', zGridPts=None)
        # interpolate
        ATx = self.interpGrid(self.alongTrackUnitVectorXInterpolator, x, y, z,
                              maskNoData=maskNoData)
        ATy = self.interpGrid(self.alongTrackUnitVectorYInterpolator, x, y, z,
                              maskNoData=maskNoData)
        if save:
            self.alongTrackUnitVectorX = ATx.astype(np.float32)
            self.alongTrackUnitVectorY = ATy.astype(np.float32)
        return ATx, ATy

    def losUnitVectorCube(self, x, y, z, maskNoData=None, save=False):
        '''
        Interpolate losUnitVector data cube

        Parameters
        ----------
        x or slantTrange : nparray
            x (G) or slant (R) coordinates.
        y : nparray
            y (G) or slant (ZD) coordinates.
        z : nparray
            Elevation coordinates.
        save: bool, optional
            If true, save the result as self.result (e.g., self.losUnitVectorX,
            self.losUnitVectorY). The default is False.

        Returns
        -------
        LOSx, LOSy : nparrays
             x and y components of the los unit vector

        '''
        # Setup interpolators if needed.
        if not hasattr(self, 'losUnitVectorXInterpolator'):
            self.setupDataCube('losUnitVectorX', zGridPts=None)
        if not hasattr(self, 'losUnitVectorYInterpolator'):
            self.setupDataCube('losUnitVectorY', zGridPts=None)
        # interpolate
        LOSx = self.interpGrid(self.losUnitVectorXInterpolator, x, y, z,
                               maskNoData=maskNoData)
        LOSy = self.interpGrid(self.losUnitVectorYInterpolator, x, y, z,
                               maskNoData=maskNoData)
        if save:
            self.losUnitVectorX = LOSx.astype(np.float32)
            self.losUnitVectorY = LOSy.astype(np.float32)
        return LOSx, LOSy

    # Gelocaton code that should not be needed given geolocation cubes
    # Most is GrIMP specfic except what is need for
    # self.SLCSceneCenterAlongTrackSpacing

    def computeAngles(self, slantRange, zeroDopplerTime, zWGS84,
                      degrees=False):
        '''
        Compute elevation angle for slantRange, zeroDoppler, zWGS84.

        Note this elevation angle is in the zero-Doppler plane orthogonal to
        the satellite track (the dataCube elevation angle includes the squint)

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

    def getSatelliteHeight(self, zdTime):
        '''
        Compute satellite height along track

        Parameters
        ----------
        zdTime : list or np.array of points
            Zero Doppler times at which to compute the satellite height.

        Returns
        -------
        hSat : list or np.array of points
            Heights corresponding to zero Doppler times.

        '''
        if not hasattr(self, 'ECEFtoLL'):
            self.ECEFtoLL = pyproj.Transformer.from_crs("EPSG:4978",
                                                        "EPSG:4326").transform
        # Get the ECEF coords
        x = self.orbit.xsv(zdTime)
        y = self.orbit.ysv(zdTime)
        z = self.orbit.zsv(zdTime)
        # Convert to geodetic lat, lon, height above ellipsoid
        latSat, lonSat, hSat = self.ECEFtoLL(x, y, z)
        return hSat

    def RTtoLatLon(self, slantRange, zdTime, z):
        '''
        Convert R/D coords to lat lon at a point. In most cases shouldn't be
        needed because this can be acomplished faster with the data cube.

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

    def RTtoECEF(self, slantRange, zdTime, z):
        '''
        Convert R/D coords to lat lon at a point. In most cases shouldn't be
        needed because this can be acomplished faster with the data cube.

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
        return self.SMLocateZDECEF(position, velocity, slantRange, z)

    def SMLocateZDECEF(self, position, velocity, slantRange, z):
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
        x, y, z
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
        # Function to optimize
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
        # return location
        return t

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
             coordinates

        '''
        t = self.SMLocateZDECEF(position, velocity, slantRange, z)
        return self.ECEFtoLL(*list(t))

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
