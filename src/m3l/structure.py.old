#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:07:12 2022

@author: flaviano.fernandes@ifpr.edu.br
"""

import logging
from m3l.utils import Constant

logging.basicConfig(
    level=logging.WARNING,
    filename="structure.log",
    format="%(asctime)s:%(levelname)s:%(message)s"
    )

class Atoms(Constant):

    atoms = list()

    def __init__(self) -> None:
    
        super().__init__()

    def filterAttr(self, attr1, attr2, value):

        var = list()
        for atom in self.atoms:
            if atom[attr2] == value:
                var.append(atom[attr1])

        return var

class Lattice(Atoms):

    def __init__(self, a, b, c):

        super().__init__()

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

