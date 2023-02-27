from math import exp, sqrt, cos, pi

class Symm2D():

    def __init__(self, atoms, prm, rs):
        self.atoms = atoms
        self.eta_prm = prm
        self.rs_prm = rs

    def gFunction(self):
        
        for atm in self.atoms:
            sum = 0.e0
            for atom in self.atoms:
                if atom['id'] > atm['id']:
                    dr = sqrt((atm['x']-atom['x'])**2+(atm['y']-atom['y'])**2+(atm['z']-atom['z'])**2)
                    sum += exp(-self.eta_prm*(dr-self.rs_prm)**2)*self.setSF(dr)
            atm['s2'] = sum

        return self.atoms

    def setSF(self, dr):

        if dr <= self.rs_prm:
            func = 0.5*(cos(pi*dr/self.rs_prm)+1.e0)
        else:
            func = 0.e0

        return func
