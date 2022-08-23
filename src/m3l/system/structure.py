#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:07:12 2022

@author: flaviano
"""

class Cell():
    def __init__(self, a, b, c):
        self.acell = a
        self.bcell = b
        self.ccell = c
        self.a1 = a
        self.a2 = 0.e0
        self.a3 = 0.e0
        self.b1 = 0.e0
        self.b2 = b
        self.b3 = 0.e0
        self.c1 = 0.e0
        self.c2 = 0.e0
        self.c3 = c

    def __str__(self):
        return ('\nCélula unitária ortorrômbica ({}, {}, {})'
                .format(self.acell, self.bcell, self.ccell))