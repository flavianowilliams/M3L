import numpy as np
from m3l import libs
from m3l.utils import Constants

class ForceField(Constants):

    params = np.array([])

    def potential(self, *args):

        self.params = np.array([args[0]])

        for i in range(1, len(args)):

            self.params = np.append(self.params, [args[i]], axis = 0)

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

    def site(charge, eps, sigma):

        return np.array([charge, eps, sigma])

