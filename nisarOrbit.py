#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 09:39:32 2024

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
from xml.dom.minidom import parseString




class nisarOrbit():
    '''
    Class to read and interpret NISAR orbit data.
    '''

    __metaclass__ = ABCMeta

    def __init__(self, h5OrbitGroup=None, XMLOrbit=None):
        '''
        Init nisarOrbit. Load orbit data if provided.
        Parameters
        ----------
        h5Orbit : h5 orbit group, optional
            orbit group from an h5 file. The default is None.
        xmlOrbit : str, optional
            XML orbit file. The default is None.

        Returns
        -------
        None.

        '''
        self.parsers = {'utc': self._parseTime, 'tai': self._parseTime,
                        'gps': self._parseTime, 'x':  self._parseFloat,
                        'y': self._parseFloat, 'z': self._parseFloat,
                        'vx':  self._parseFloat, 'vy':  self._parseFloat,
                        'vz':  self._parseFloat}
        if h5OrbitGroup is not None and XMLOrbit is not None:
            self.printError('Cannot specify both h5Orbit and xmlOrbit')
            return
        print(h5OrbitGroup, XMLOrbit)
        if h5OrbitGroup is not None:
            self.parseStateVectorsH5(h5OrbitGroup)
        if XMLOrbit is not None:
            self.parseStateVectorsXML(XMLOrbit)
        #
        # format
        self.formatStateVectors()

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

    def parseStateVectorsH5(self, h5OrbitGroup):
        '''
        Parse the state vectors. Store with names compatiable with geodat
        files.
        Parameters
        ----------
        h5Orbit : h5 orbit group
            Orbit group from an h5 file.
        Returns
        -------
        None.

        '''
        # Get time array
        self.time = np.array(h5OrbitGroup['time'])
        #
        # use it to compute other time related values
        self.NumberOfStateVectors = len(self.time)
        self.TimeOfFirstStateVector = self.time[0]
        self.StateVectorInterval = self.time[1] - self.time[0]
        #
        # Get the orbit type
        self.stateVectorType = np.array(h5OrbitGroup['orbitType']
                                        ).item().decode()
        #
        # Get position and velocity
        self.position = np.array(h5OrbitGroup['position'])
        self.velocity = np.array(h5OrbitGroup['velocity'])
        #
        # Create interpolationrs pos and vel for state vectors
        for i, pos, vel in zip(range(0, 3),
                               ['xsv', 'ysv', 'zsv'],
                               ['vxsv', 'vysv', 'vzsv']):
            setattr(self, pos,
                    scipy.interpolate.RegularGridInterpolator(
                        [self.time], self.position[:, i], method='quintic'))
            setattr(self, vel,
                    scipy.interpolate.RegularGridInterpolator(
                        [self.time], self.velocity[:, i], method='quintic'))

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

    def _parseFloat(self, textString):
        return np.float64(textString)

    def _parseTime(self, textString):
        textString = textString[0:min(26, len(textString))]
        #
        return datetime.strptime(textString, '%Y-%m-%dT%H:%M:%S.%f')

    def _checkChildren(self, xmlData, level, sv):
        if not xmlData.hasChildNodes():
            if xmlData.nodeType == xmlData.TEXT_NODE:
                if level == 4:
                    if xmlData.parentNode.tagName in self.parsers.keys():
                        sv[xmlData.parentNode.tagName].append(
                            self.parsers[xmlData.parentNode.tagName](
                                xmlData.data)
                            )

        for children in xmlData.childNodes:
            self._checkChildren(children, level+1, sv)

    def parseStateVectorsXML(self, XMLFile, referenceDate=None):
        '''
        Parse the state vector data from the xml file

        Parameters
        ----------
        XMLFile : str
            Name of xml file with state vector data.

        Returns
        -------
        sv : dict
            dictionary wit time (tai, utc, gps), position(x,y,z), and
            velocity (vx, vy, vy) data.
        '''
        #
        # Initial state vectors
        sv = {'tai': [], 'utc': [], 'gps': [],
              'x': [], 'y': [], 'z': [], 'vx': [], 'vy': [], 'vz': []}
        #
        # Open the xml file
        with open(XMLFile, 'r') as fp:
            for line in fp:
                # Avoid any blank lines
                if len(line) > 200:
                    break
        #
        # Loop to find sv data, then quit loop
        for xmlData in parseString(line).childNodes:
            self._checkChildren(xmlData, 0, sv)
            if len(sv['x']) > 5:
                break
        #
        if referenceDate is None:
            referenceDate = sv['utc'][0].replace(hour=0, minute=0,
                                                 second=0, microsecond=0)
        #
        print('first', sv['utc'][0])
        self.NumberOfStateVectors = len(sv['utc'])
        self.TimeOfFirstStateVector = (sv['utc'][0] -
                                       referenceDate).total_seconds()
        self.StateVectorInterval = (sv['utc'][1] -
                                    sv['utc'][0]).total_seconds()
        print(self.TimeOfFirstStateVector , self.StateVectorInterval)
        #
        # Get the orbit type
        #self.stateVectorType = np.array(h5OrbitGroup['orbitType']
        #                                 ).item().decode()
        #
        # Get position and velocity
        self.position = np.column_stack([sv['x'], sv['y'], sv['z']])
        self.velocity =np.column_stack([sv['vx'], sv['vy'], sv['vz']])
