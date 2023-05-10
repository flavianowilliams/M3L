from math import exp, sqrt, cos, pi

class Symm2D():

    def __init__(self, atoms, atm, eta_prm, rs_prm):

        self.atoms = atoms
        self.eta_prm = eta_prm
        self.rs_prm = rs_prm
        self.atm = atm

        self.gFunction()

    def gFunction(self):
        
        sum = 0.e0
        for atom in self.atoms:
            if atom != self.atm:
                dr = sqrt((self.atm['x']-atom['x'])**2+(self.atm['y']-atom['y'])**2+(self.atm['z']-atom['z'])**2)
                sum += exp(-self.eta_prm*(dr-self.rs_prm)**2)*self.setSF(dr)
        self.atm['s2'] = sum

    def setSF(self, dr):

        if dr <= self.rs_prm:
            func = 0.5*(cos(pi*dr/self.rs_prm)+1.e0)
        else:
            func = 0.e0

        return func
