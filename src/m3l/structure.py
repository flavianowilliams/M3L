#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:07:12 2022

@author: flaviano
"""

import logging
from m3l.utils import Constant

logging.basicConfig(
    level=logging.WARNING,
    filename="structure.log",
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

class Atoms(Constant):

    def __init__(self):

        self.atoms = list()

class Lattice(Atoms):

    def __init__(self, a, b, c):

        self.__acell = a
        self.__bcell = b
        self.__ccell = c

    def setAcell(self, a):
        self.__acell = a

    def getAcell(self):
        return self.__acell

    def setBcell(self, b):
        self.__bcell = b

    def getBcell(self):
        return self.__bcell

    def setCcell(self, c):
        self.__ccell = c

    def getCcell(self):
        return self.__ccell

    def setVolume(self):
        self.__volume = {'value': self.__acell*self.__bcell*self.__ccell, 'unit': 'A³'}

    def getVolume(self):      
        return self.__volume

