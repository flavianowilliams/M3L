from m3l.structure import Lattice
from m3l.symmetry import Symm2D
import logging
import sys

class System(Lattice):

    def __init__(self, a, b, c, filename, eta_prm, rs_prm):

        self.setSym()
        self.eta_prm = eta_prm
        self.rs_prm = rs_prm

    def setXYZ(self, filename):
        with open(filename, 'r') as xyz_file:
            self.__natoms = xyz_file.readline()
            self.__natoms = int(self.__natoms)
            xyz_file.readline()
            self.atoms = list()
            for i in range(self.__natoms):
                p1, p2, p3, p4 = xyz_file.readline().split(maxsplit=3)
                p2 = float(p2)
                p3 = float(p3)
                p4 = float(p4)
                self.atoms.append({'id': i, 'atom': p1, 'x': p2, 'y': p3, 'z': p4, 'energy': 0.e0, 's2': 0.e0})

    def setSym(self):
        for atom in self.atoms:
            atom = Symm2D(self.atoms, atom, self.eta_prm, self.rs_prm)

    def setCCP(self):
        if self.getAcell() > 0.0 and self.getBcell() > 0.0 and self.getCcell() > 0.0:
            for at in self.atoms:
                at['x'] = at['x']+self.getAcell()*int(at['x']/self.getAcell())
                at['y'] = at['y']+self.getBcell()*int(at['y']/self.getBcell())
                at['z'] = at['z']+self.getCcell()*int(at['z']/self.getCcell())
        else:
            logging.critical('Constant lattice must be larger than zero!')
            sys.exit()

    def getNatoms(self):
        return self.__natoms
