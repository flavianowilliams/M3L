#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:07:12 2022

@author: flaviano
"""

import numpy as np

class Lattice():
    def __init__(self, a, b, c):
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
        self.__volume = self.__acell*self.__bcell*self.__ccell

    def getVolume(self):
        return self.__volume

class System(Lattice):
    def __init__(self, a, b, c, filename):
        self.setAcell(a)
        self.setBcell(b)
        self.setCcell(c)
        self.setVolume()
        self.setXYZ(filename)

    def setXYZ(self, filename):
        with open(filename, 'r') as xyz_file:
            self.__natoms = xyz_file.readline()
            self.__natoms = int(self.__natoms)
            xyz_file.readline()
            self.atoms = list()
            for i in range(self.__natoms):
                p1, p2, p3, p4 = xyz_file.readline().split(maxsplit=3)
                p2 = float(p2)
                p3 = float(p3)
                p4 = float(p4)
                self.atoms.append({'id': i, 'atom': p1, 'x': p2, 'y': p3, 'z': p4})

    def getNatoms(self):
        return self.__natoms

    def reciprocal_lattice(self):
        for at in self.atoms:
            at['x'] = at['x']+at['x']*int(at['x']/self.getAcell())
            at['y'] = at['y']+at['y']*int(at['y']/self.getBcell())
            at['z'] = at['z']+at['z']*int(at['z']/self.getCcell())

    def __str__(self):
        return (
            '\nCélula unitária ortorrômbica\n\n'
            +'Volume: {} A**3\n\n'.format(self.getVolume())
            +'Total: {} átomos\n'.format(self.getNatoms())
        )
