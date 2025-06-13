import numpy as np
from m3l.utils import Constants

# classe responsavel pela definicao do campo de força 

class ForceField2(Constants):

    van_der_waals = []
    electrostatic = None

    def intermolecular(self, van_der_waals, electrostatic, rvdw = 0.0e0, rcoul = 0.0e0):

        self.rvdw = rvdw
        self.rcoul = rcoul
        self.van_der_waals = van_der_waals
        self.electrostatic = electrostatic

    def hook(self):

        vdw = self.van_der_waals
        electrostatic = self.electrostatic

        self.intermolecular(vdw, electrostatic)

        self.nvdw = len(self.van_der_waals)

    def hook_output(self):

        return (self.nvdw, self.rvdw, self.rcoul, self.van_der_waals)

    def __call__(self):
        
        self.hook()

        return self.hook_output()

    def convertUnits(self):

        for item in self.van_der_waals:
            
            item[2] = item[2]*self.ECONV
            item[2] = item[2].item()
            item[3] = item[3]*self.ACONV
            item[3] = item[3].item()

# classe responsavel por definir as interaçoes de Van der Waals 

class Intermolecular():

    def site(s1, s2, eps, sigma):

        list_ = [s1, s2, eps, sigma]

        return list_

    def sites(*args):

        list_ = []
        for item in args:
            list_.append(item)

        return list_

# classe responsavel para definicao das moleculas

class Molecules(Constants):

    def molecule(self, nsite):

        return np.array([nsite], dtype = np.float64)

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

class Intramolecular():

    def bond(n1, n2, k, r0):

        utls = Constants()
        k = k*utls.ECONV/(utls.ACONV)**2 

        return np.array([n1 , n2, k, r0], dtype = np.float64)
