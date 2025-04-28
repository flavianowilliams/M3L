import numpy as np
from m3l import libs
from m3l.utils import Constants

class ForceField(Constants):

    prms = np.array([np.zeros(4)], dtype = np.float64)
    vdw = np.array([])

    def structure(self, *args):

        for arg in args:

            if self.vdw.size == 0:

                self.vdw = np.array(arg, dtype = np.float64)

            else:

                self.vdw = np.append(self.vdw, arg, axis = 0)

    def intermolecular(self, **kwargs):

        self.prms[0, 0] = kwargs['rvdw']

    def molecule(self, *args):

        return args

    def interaction(self, system, potential):

        libs.rx = system.rx
        libs.ry = system.ry
        libs.rz = system.rz
        libs.natom = len(system.zat)
        libs.cell = system.cell
        libs.params = potential
        libs.forces()

        system.epotential = libs.energy

        self.fx = libs.fx 
        self.fy = libs.fy 
        self.fz = libs.fz 
        self.ea = libs.ea 

        return system

    def hook(self):

        self.params = np.append(self.prms, self.vdw, axis = 0)

    def __call__(self):

        self.hook()

        return self.params

class Intermolecular(Constants):

    def site(indx, charge, eps, sigma):

        return np.array([indx, charge, eps, sigma], dtype = np.float64)

