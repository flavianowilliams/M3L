#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:07:12 2022

@author: flaviano
"""

import logging
import sys
from m3l.utils import Constant
from m3l.symmetry import Symm2D

logging.basicConfig(
    level=logging.WARNING,
    filename="structure.log",
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

class Lattice(Constant):
    def __init__(self, a, b, c):

        super().__init__()

        self.__acell = a
        self.__bcell = b
        self.__ccell = c
        self.__a1 = a
        self.__a2 = 0.e0
        self.__a3 = 0.e0
        self.__b1 = 0.e0
        self.__b2 = b
        self.__b3 = 0.e0
        self.__c1 = 0.e0
        self.__c2 = 0.e0
        self.__c3 = c

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
