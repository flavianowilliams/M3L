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

class System(Lattice):

    def __init__(self, a, b, c, filename, eta_prm, rs_prm):

        self.setAcell(a)
        self.setBcell(b)
        self.setCcell(c)
        self.eta_prm = eta_prm
        self.rs_prm = rs_prm
        self.setXYZ(filename)
        self.convertUnits()
        self.setCCP()
        self.setSym()
        self.setVolume()

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
                self.atoms.append({'id': i, 'atom': p1, 'x': p2, 'y': p3, 'z': p4, 'energy': 0.e0, 's2': 0.e0})

    def setSym(self):
        for atom in self.atoms:
            atom = Symm2D(self.atoms, atom, self.eta_prm, self.rs_prm)

    def setCCP(self):
        if self.getAcell() > 0.0 and self.getBcell() > 0.0 and self.getCcell() > 0.0:
            for at in self.atoms:
                at['x'] = at['x']+self.getAcell()*int(at['x']/self.getAcell())
                at['y'] = at['y']+self.getBcell()*int(at['y']/self.getBcell())
                at['z'] = at['z']+self.getCcell()*int(at['z']/self.getCcell())
        else:
            logging.critical('Constant lattice must be larger than zero!')
            sys.quit()

    def getNatoms(self):
        return self.__natoms
