#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 21:05:45 2022

@author: flaviano
"""

from pathlib import Path
import sys

path_root = Path(__file__).parents[0]
sys.path=list()
sys.path.insert(0, str(path_root))

import settings
from system import structure

def teste():
    return print(structure.mensagem)