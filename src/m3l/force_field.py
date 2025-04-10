import numpy as np
from m3l import libs
from m3l.utils import Constants

class ForceField(Constants):

    params = np.array([], dtype = np.float64)

    def potential(self, *args):

        list_ = []
        for arg in args:

            list_.append(arg[0])

        list_.sort()

        for i, arg in enumerate(args):

            indx = list_.index(arg[0])

            if self.params.size == 0:

                self.params = np.array([args[indx]], dtype = np.float64)
            
            else:

                self.params = np.append(self.params, [args[indx]], axis = 0)

#        self.params = np.delete(self.params, 0, axis = 1)

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

    def __call__(self):

        return self.params

class Intermolecular(Constants):

    def site(indx, charge, eps, sigma):

        return np.array([indx, charge, eps, sigma], dtype = np.float64)

