import numpy as np
from m3l.utils import Constants

class ForceField(Constants):

    prms = np.array([np.zeros(4)], dtype = np.float64)
    vdw = np.array([])
    bond = np.array([])

    def structure(self, *args):

        for arg in args:

            if self.vdw.size == 0:

                self.vdw = np.array(arg, dtype = np.float64)

            else:

                self.vdw = np.append(self.vdw, arg, axis = 0)

    def intermolecular(self, **kwargs):

        self.prms[0, 0] = kwargs['rvdw']*self.ACONV
        self.prms[0, 1] = kwargs['rcoul']*self.ACONV

    def molecule(self, *args):

        return args

    def hook(self):

        self.params = self.prms

        if self.vdw.size > 0:

            self.params = np.append(self.params, self.vdw, axis = 0)

    def __call__(self):

        self.hook()

        return self.params

class Intermolecular():

    def site(indx, charge, eps, sigma):

        utls = Constants()
        eps = eps*utls.ECONV
        sigma = sigma*utls.ACONV

        return np.array([indx, charge, eps, sigma], dtype = np.float64)

class Intramolecular():

    def bond(n1, n2, k, r0):

        utls = Constants()
        k = k*utls.ECONV/(utls.ACONV)**2 

        return np.array([n1 , n2, k, r0], dtype = np.float64)
