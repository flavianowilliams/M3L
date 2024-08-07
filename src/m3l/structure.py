import numpy as np
import json
from m3l.utils import Conversion
from numpy._core.multiarray import dtype

class Atom(Conversion):

    atoms = []

    def __init__(self, zat):
        super().__init__()

        self.zat = zat

    def setMass(self):

        if int(self.zat) == 1:
            self.mass = 1.008
        elif int(self.zat) == 18:
            self.mass = 39.948

        return self.mass/self.MCONV

class System():

    def __init__(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.cell = np.array(json_file['cell'], dtype=np.float32)
            self.atoms = np.array(json_file['atoms'], dtype=np.float32)
            self.temperature = np.array(json_file['thermodynamic'][0], dtype=np.float32)
            self.pressure = np.array(json_file['thermodynamic'][1], dtype=np.float32)

    def setNAtoms(self):
        self.natoms = len(self.atoms)
        return self.natoms

class Symmetry(System):

    def __init__(self):

        self.eta_prm = 0.0
        self.rs_prm = 0.0
        self.rc_prm = 6.0

    def symm2D(self):
        
        for atm in self.atoms:
            sum_ = 0.e0
            for atom in self.atoms:
                if atom[0] != atm[0]:
                    dr = np.sqrt((atm[4]-atom[4])**2+(atm[5]-atom[5])**2+(atm[6]-atom[6])**2)
                    sum_ += np.exp(-self.eta_prm*(dr-self.rs_prm)**2)*self.setSF(dr)
            return sum_

    def setSF(self, dr):

        if dr <= self.rc_prm:
            func = 0.5*(np.cos(np.pi*dr/self.rc_prm)+1.e0)
        else:
            func = 0.e0

        return func

