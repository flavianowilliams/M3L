from math import exp, sqrt, cos, pi
from m3l.utils import Constant

class Symm2D(Constant):

    def __init__(self, atoms, prm, rs):
    
        self.sym_coord = list()
        self.eta_prm = prm
        self.rs_prm = rs

        self.gFunction(atoms)

    def gFunction(self, atoms):
        
        for atm in atoms:
            sum = 0.e0
            for atom in atoms:
                if atom['id'] > atm['id']:
                    dr = sqrt((atm['x']-atom['x'])**2+(atm['y']-atom['y'])**2+(atm['z']-atom['z'])**2)
                    sum += exp(-self.eta_prm*(dr-self.rs_prm)**2)*self.setSF(dr)
            atm['s2'] = sum
            self.sym_coord.append(atm)

        return self.sym_coord

    def setSF(self, dr):

        if dr <= self.rs_prm:
            func = 0.5*(cos(pi*dr/self.rs_prm)+1.e0)
        else:
            func = 0.e0

        return func
