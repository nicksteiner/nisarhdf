#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 08:36:38 2024

@author: ian
"""
from nisarHDF import nisarBaseHDF
from osgeo import gdal
import os
import numpy as np
import rioxarray


class nisarROFFHDF(nisarBaseHDF):
    '''
    This class creates objects to work with nisar RUNWHDF images.
    '''

    def __init__(self,  sar='LSAR', frequency='frequencyA', polarization='HH',
                 layer='layer1', byteOrder='LSB'):
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
        layer : str, optional
            Layer to use. The default is 'layer1'.

        Returns
        -------
        None.

        '''
        nisarBaseHDF.__init__(self,
                              sar=sar,
                              product='ROFF',
                              frequency=frequency,
                              productType='pixelOffsets',
                              polarization=polarization,
                              layer=layer,
                              productData='alongTrackOffset',
                              bands='swaths',
                              bytOrder=byteOrder)

    def parseParams(self, secondGeodat=None, **kwds):
        '''
        Parse all the params needed for offsets

        Returns
        -------
        None.

        '''
        self.isSecondary = False
        self.getOffsetParams()
        # self.getNumberOfLooks()
        self.getSize(offsets=True)
        self.getSingleLookPixelSizeOffsets()
        self.getEPSG()

    def writeOffsetsVrt(self, vrtFile, sourceFiles, geodat1, geodat2, mask,
                        Image1=None, Image2=None, scaleFactor=1,
                        descriptions=None, byteOrder='LSB'):
        metaData = {}
        #
        if descriptions is None:
            descriptions = [os.path.basename(x) for x in sourceFiles]
        #
        for key in ['a0', 'r0', 'deltaA', 'deltaR']:
            metaData[key] = getattr(self, key)

        for key, var in zip(['geo1', 'geo2', 'mask', 'Image1', 'Image2'],
                            [geodat1, geodat2, mask, Image1, Image2]):
            if var is not None:
                metaData[key] = var
        #
        self._writeVrt(vrtFile,
                       sourceFiles,
                       descriptions,
                       byteOrder=byteOrder,
                       eType=gdal.GDT_Float32,
                       geoTransform=[-0.5, 1., 0., -0.5, 0., 1.],
                       metaData=metaData,
                       setSRS=False,
                       noDataValue=-2.e9


                       )

    def writeOffsetsDatFile(self, datFile, geodat1=None, geodat2=None):
        '''
        Write an offsets data file

        Parameters
        ----------
        datFile : TYPE
            DESCRIPTION.
        geodat1 : TYPE
            DESCRIPTION.
        geodat2 : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        with open(datFile, 'w') as fp:
            print(self.a0, self.r0,
                  self.MLRangeSize, self.MLAzimuthSize,
                  self.deltaR, self.deltaA,
                  file=fp)
            if geodat1 is None:
                return
            if geodat2 is None:
                geodat2 = geodat1.replace('.in', '.secondary'
                                          ).replace('.geojson', '.secondary')
            print(geodat1, geodat2, file=fp)

    def getLayers(self, layerTypes, layers=['layer1', 'layer2', 'layer3']):
        '''
        Get Layers for a given layerNames

        Parameters
        ----------
        layerName : TYPE
            DESCRIPTION.
        layers : TYPE, optional
            DESCRIPTION. The default is ['layer1', 'layer2', 'layer3'].

        Returns
        -------
        None.

        '''
        h5Data = \
            self.h5[self.product]['swaths'][self.frequency][self.productType]
        for layerType in layerTypes:
            layerData = []
            for layer in layers:
                layerData.append(h5Data[self.polarization][layer][layerType])
            setattr(self, layerType, np.stack(layerData))

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
        # Load the product fields if not already done
        self.getProducts([filterField, 'alongTrackOffset', 'slantRangeOffset'])
        # Loop over layers
        qaData = getattr(self, filterField)
        azOffsets = getattr(self, 'alongTrackOffset')
        rgOffsets = getattr(self, 'slantRangeOffset')
        for layer in layers:
            qaLayer = qaData[layer-1]
            bad = qaLayer < thresholds[layer-1]
            azOffsets[layer-1][bad] = np.nan
            rgOffsets[layer-1][bad] = np.nan

    def applyMask(self, maskFile, layers=[3], maskValue=1):
        '''
        Apply mask contained in file

        Parameters
        ----------
        maskFile : str
            vrt file with mask (offsets.mask.vrt)
        t
        layers : list, optional
            List of layers to process. The default is [3].

        Returns
        -------
        None.

        '''
        #
        mask = rioxarray.open_rasterio(maskFile)
        masked = np.squeeze(mask.data) > maskValue
        # Load the product fields if not already done
        self.getProducts(['alongTrackOffset', 'slantRangeOffset'])
        # Loop over layers
        for layer in layers:
            self.alongTrackOffset[layer-1][masked] = np.nan
            self.slantRangeOffset[layer-1][masked] = np.nan

    def getProducts(self, productFields):
        '''
         Load the product specificed product fields if not already done.

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

    def writeOffsets(self, filenameRoot,
                     geojsonName=None,
                     geojsonNameSecondary=None,
                     productFields=['slantRangeOffset',
                                    'alongTrackOffset',
                                    'correlationSurfacePeak'],
                     suffixes=['.dr', '.da', '.cc'],
                     dataTypes=['>f4', '>f4', '>f4'],
                     layers=[1, 2, 3], matchTypes=[1, 2, 3]):
        '''
        Breakout offsets and write as layers.

        Parameters
        ----------
        filenameRoot : str
            Basename for offsets (e.g., offsets.dr, .da...).
        geojsonName : str optional
            geojson for redference image. The default is None.
        geojsonNameSecondary : TYPE, optional
            geojson for secondary. The default is None.
        productFields : list of str, optional
            Product Fields to save. The default is ['alongTrackOffset',
                                                    'slantRangeOffset','snr'].
        suffixes : list of str, optional
            suffixes for output fields. Should correspond to product fields.
            The default is ['.da', '.dr', '.snr'].
        dataTypes : str, optional
            Format for output. The default is byteswap ['>f4', '>f4', '>f4'].
        layers : list of ints, optional
            Layers to save. The default is [1, 2, 3].

        Returns
        -------
        None.

        '''
        mt = None
        # Make sure the data have been read
        self.getProducts(productFields)
        # Scale factors for conversion to pixels
        scaleFactors = {'slantRangeOffset': 1./self.SLCRangePixelSize,
                        'alongTrackOffset': 1./self.SLCAzimuthPixelSize}
        #
        # Look over layers
        for layer, matchType in zip(layers, matchTypes):
            sourceFiles = []
            descriptions = []
            for productField, suffix, dataType in \
                    zip(productFields, suffixes, dataTypes):
                data = getattr(self, productField)
                sourceFile = f'{filenameRoot}.layer{layer}{suffix}'
                sourceFiles.append(os.path.basename(sourceFile))
                descriptions.append(f'{productField} {layer}')
                if productField in scaleFactors:
                    scaleFactor = scaleFactors[productField]
                else:
                    scaleFactor = 1.
                self._writeImageData(sourceFile,
                                     np.squeeze(data[layer-1, :, :]
                                                ) * scaleFactor, dataType)
                #
                # If along track offsets save match type
                if productField == 'alongTrackOffset':
                    mt = np.zeros(np.squeeze(data[layer-1, :, :]).shape,
                                  'byte')
                    mt[np.squeeze(data[layer-1, :, :]) > -2.e8] = matchType
                    self._writeImageData(f'{filenameRoot}.layer{layer}.mt',
                                         mt, 'u1')

            # Vrt file for the layer
            vrtFile = f'{filenameRoot}.layer{layer}.vrt'
            byteOrder = 'LSB'
            if '>' in dataType:
                byteOrder = 'MSB'
            # Save the vrt
            self.writeOffsetsVrt(vrtFile, sourceFiles, geojsonName,
                                 geojsonNameSecondary, None,
                                 Image1=None,
                                 Image2=None,
                                 scaleFactor=1,
                                 descriptions=descriptions,
                                 byteOrder=byteOrder)
            #
            self.writeOffsetsDatFile(f'{filenameRoot}.layer{layer}.dat',
                                     geodat1=geojsonName,
                                     geodat2=geojsonNameSecondary)
            descriptions = [f'matchType {layer}']
            #
            if mt is not None:
                self._writeVrt(f'{filenameRoot}.layer{layer}.mt.vrt',
                               [f'{filenameRoot}.layer{layer}.mt'],
                               descriptions,
                               eType=gdal.GDT_Byte,
                               geoTransform=[-0.5, 1., 0., -0.5, 0., 1.],
                               setSRS=False,
                               noDataValue=0)

    def writeData(self, filenameRoot, productField,  layers=[1, 2, 3],
                  geojsonName=None, geojsonNameSecondary=None,
                  dataType='>f4', metaData=None, datFile=None,
                  suffix=''):
        '''
       Write offset data
       Parameters
       ----------
       filename : TYPE
           DESCRIPTION.
       productField : TYPE
           DESCRIPTION.
       layers : TYPE, optional
           DESCRIPTION. The default is ['layer1', 'layer2', 'layer3'].
       geojsonName : str, optional
           DESCRIPTION. The default is None.
       geojsonNameSecondary : TYPE, optional
           DESCRIPTION. The default is None.
       dataType : str, optional
           DESCRIPTION. The default is '>f4'.
       metaData : TYPE, optional
           DESCRIPTION. The default is None.
       suffix : str, optional
           DESCRIPTION. The default is None.

       Returns
       -------
       None.

       '''
        #
        # Write the image data to a floating point file
        if not hasattr(self, productField):
            print(f'getting layers for {productField}')
            self.getLayers([productField])

        data = getattr(self, productField)
        sourceFiles = []
        descriptions = []
        for layer in layers:
            sourceFile = f'{filenameRoot}.layer{layer}{suffix}'
            print(f'writing layer {sourceFile}')
            print(np.min(data), np.max(data))
            self._writeImageData(sourceFile,
                                 np.squeeze(data[layer-1, :, :]), dataType)
            sourceFiles.append(os.path.basename(sourceFile))
            descriptions.append(f'{productField} {layer}')
        vrtFile = f'{filenameRoot}{suffix}.vrt'
        #
        byteOrder = 'LSB'
        if '>' in dataType:
            byteOrder = 'MSB'
        #
        self.writeOffsetsVrt(vrtFile, sourceFiles,
                             geojsonName, geojsonNameSecondary, None,
                             Image1=None, Image2=None,
                             scaleFactor=1, descriptions=descriptions,
                             byteOrder=byteOrder)
        #
        if datFile is not None:
            self.writeOffsetsDatFile(datFile,
                                     geodat1=geojsonName,
                                     geodat2=geojsonNameSecondary)
        #
        # write geodat in same dir as filename with geodatNRxNA unless other
        # if includeGeojson:
        #     self.writeGeodatGeojson(filename=geojsonName,
        #                             path=os.path.dirname(filename),
        #                             secondary=secondary)
        # #
        # # Write an accompanyting vrt file if requested filename.vrt
        # if includeVRT:
        #     if '>' in dataType:
        #         byteOrder = 'MSB'
        #     else:
        #         byteOrder = 'LSB'
        #     # pixel centered, lower left corner, pixel units
        #     geoTransform = [-0.5, 1., 0., -0.5, 0., 1.]
        #     self._writeVrt(f'{filename}.vrt',
        #                    [os.path.basename(filename)],
        #                    [productField],
        #                    metaData=metaData,
        #                    byteOrder=byteOrder,
        #                    setSRS=False,
        #                    geoTransform=geoTransform)
