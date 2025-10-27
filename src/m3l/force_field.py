import numpy as np
from m3l.utils import Constants

# classe responsavel pela definicao do campo de força 

class ForceField(Constants):

    van_der_waals = []
    electrostatic = None

    def sites(self, *args):

        list1_ = []
        list2_ = []
        for item in args:
            list1_.append(item[0])
            list2_.append(item[1])

        self.van_der_waals = [list1_, list2_]

        return self.van_der_waals

    def intermolecular(self, van_der_waals, electrostatic, rvdw = 0.0e0, rcoul = 0.0e0):

        self.van_der_waals = van_der_waals
        self.electrostatic = electrostatic
        self.rvdw = rvdw
        self.rcoul = rcoul

    def hook(self):

        vdw = self.van_der_waals
        electrostatic = self.electrostatic

        self.intermolecular(vdw, electrostatic)

        self.nvdw = len(self.van_der_waals[0])

    def hook_output(self):

        return (self.nvdw, self.rvdw, self.rcoul, self.van_der_waals)

    def __call__(self):
        
        self.hook()

        return self.hook_output()

    def convertUnits(self):

        for item in self.van_der_waals[1]:
            
            item[0] = item[0]*self.ECONV
            item[0] = item[0].item()
            item[1] = item[1]*self.ACONV
            item[1] = item[1].item()

# classe responsavel por definir as interaçoes de Van der Waals 

class Intermolecular():

    def site(s1, s2, eps, sigma):

        return [s1, s2], [eps, sigma]

# classe responsavel para definicao das moleculas

class Molecules(Constants):

    def molecule(self, nsite):

        return np.array([nsite], dtype = np.float64)

class Intramolecular():

    def bond(n1, n2, k, r0):

        utls = Constants()
        k = k*utls.ECONV/(utls.ACONV)**2 

        return np.array([n1 , n2, k, r0], dtype = np.float64)
