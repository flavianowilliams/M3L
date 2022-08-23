#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:07:12 2022

@author: flaviano
"""

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

    def setVolume(self):
        self.__volume = self.__acell*self.__bcell*self.__ccell
        return self.__volume

    def getVolume(self):
        return self.__volume

    def __str__(self):
        return ('\nCélula unitária ortorrômbica')