import numpy as np
import json
from m3l.utils import Constants

class Atom(Constants):

    def setMass(self, zatom):

        mass = 0.0

        if int(zatom) == 1:
            mass = 1.008
        elif int(zatom) == 18:
            mass = 39.948

        return mass/self.MCONV

    def setZNumber(self, atom):
        
        zatom = None

        if atom == 'H':
            zatom = 1.0 
        elif atom == 'Ar':
            zatom = 18.0

        return zatom

class System2(Constants):

    def __init__(self, cell, temperature, friction, pressure, epotential, ekinetic, atoms):
        
        self.cell = np.array(cell)
        self.temperature = np.array(temperature)
        self.friction = np.array(friction)
        self.pressure = np.array(pressure)
        self.epotential = np.array(epotential)
        self.ekinetic = np.array(ekinetic)
        self.atoms = np.array(atoms)

    def convertUnitsInv(self):

        self.epotential = self.epotential*self.ECONV
        self.ekinetic = self.ekinetic*self.ECONV
        self.temperature = self.temperature*self.TEMPCONV
        self.friction = self.friction*(self.ECONV*self.TIMECONV/self.MCONV)

        self.cell[0] = self.cell[0]*self.ACONV
        self.cell[1] = self.cell[1]*self.ACONV
        self.cell[2] = self.cell[2]*self.ACONV

        for atom in self.atoms:
            atom[1] = atom[1]*self.ACONV
            atom[2] = atom[2]*self.ACONV
            atom[3] = atom[3]*self.ACONV
            atom[4] = atom[4]*(self.ACONV/self.TIMECONV)
            atom[5] = atom[5]*(self.ACONV/self.TIMECONV)
            atom[6] = atom[6]*(self.ACONV/self.TIMECONV)
            atom[7] = atom[7]*(self.ECONV/self.ACONV)
            atom[8] = atom[8]*(self.ECONV/self.ACONV)
            atom[9] = atom[9]*(self.ECONV/self.ACONV)
            atom[10] = atom[10]*self.ECONV

    def save(self, filename = 'system.json'):

        dictfile = {
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.friction.item(),
                    self.pressure.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'atoms': self.atoms.tolist() 
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

class System(Constants):

    def loadSystem(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.cell = np.array(json_file['cell'], dtype=np.float32)
            self.atoms = np.array(json_file['atoms'], dtype=np.float32)
            self.temperature = np.array(json_file['thermodynamic'][0], dtype=np.float32)
            self.friction = np.array(json_file['thermodynamic'][1], dtype=np.float32)
            self.pressure = np.array(json_file['thermodynamic'][2], dtype=np.float32)
            self.epotential = np.array(json_file['thermodynamic'][3], dtype=np.float32)
            self.ekinetic = np.array(json_file['thermodynamic'][4], dtype=np.float32)

    def setNAtoms(self):
        self.natoms = len(self.atoms)
        return self.natoms

    def setSystem(self, temperature, pressure, cell, xyz_file):

        self.temperature = np.array(temperature)
        self.friction = np.array(0.0)
        self.pressure = np.array(pressure)
        self.cell = np.array(cell)
        self.epotential = np.array(0.0)

        with open(xyz_file, 'r') as filename:
            natoms = filename.readline()
            natoms = int(natoms)
            self.atoms = np.zeros(11*natoms).reshape(natoms, 11)
            next(filename)
            for i in range(natoms):
                ats, x, y, z = filename.readline().split(maxsplit=3)
                self.atoms[i][0] = Atom().setZNumber(ats)
                self.atoms[i][1] = float(x)
                self.atoms[i][2] = float(y)
                self.atoms[i][3] = float(z)

        self.setVelocity()

        for atom in self.atoms:
            atom[10] = 0.0

    def setVelocity(self):

        nfree = 3*(len(self.atoms)-1)

        for atom in self.atoms:
            mass = Atom().setMass(atom[0])
            atom[4] = np.sqrt(nfree*self.temperature/(6.0*mass))
            atom[5] = np.sqrt(nfree*self.temperature/(6.0*mass)) 
            atom[6] = np.sqrt(nfree*self.temperature/(6.0*mass))

        self.ekinetic = np.array(1.5*nfree*self.KB*self.temperature)

    def convertUnits(self):

        self.epotential = self.epotential/self.ECONV
        self.ekinetic = self.ekinetic/self.ECONV
        self.temperature = self.temperature/self.TEMPCONV
        self.friction = self.friction/(self.ECONV*self.TIMECONV/self.MCONV)

        self.cell[0] = self.cell[0]/self.ACONV
        self.cell[1] = self.cell[1]/self.ACONV
        self.cell[2] = self.cell[2]/self.ACONV

        for atom in self.atoms:
            atom[1] = atom[1]/self.ACONV
            atom[2] = atom[2]/self.ACONV
            atom[3] = atom[3]/self.ACONV
            atom[4] = atom[4]/(self.ACONV/self.TIMECONV)
            atom[5] = atom[5]/(self.ACONV/self.TIMECONV)
            atom[6] = atom[6]/(self.ACONV/self.TIMECONV)
            atom[7] = atom[7]/(self.ECONV/self.ACONV)
            atom[8] = atom[8]/(self.ECONV/self.ACONV)
            atom[9] = atom[9]/(self.ECONV/self.ACONV)
            atom[10] = atom[10]/self.ECONV

    def convertUnitsInv(self):

        self.epotential = self.epotential*self.ECONV
        self.ekinetic = self.ekinetic*self.ECONV
        self.temperature = self.temperature*self.TEMPCONV
        self.friction = self.friction*(self.ECONV*self.TIMECONV/self.MCONV)

        self.cell[0] = self.cell[0]*self.ACONV
        self.cell[1] = self.cell[1]*self.ACONV
        self.cell[2] = self.cell[2]*self.ACONV

        for atom in self.atoms:
            atom[1] = atom[1]*self.ACONV
            atom[2] = atom[2]*self.ACONV
            atom[3] = atom[3]*self.ACONV
            atom[4] = atom[4]*(self.ACONV/self.TIMECONV)
            atom[5] = atom[5]*(self.ACONV/self.TIMECONV)
            atom[6] = atom[6]*(self.ACONV/self.TIMECONV)
            atom[7] = atom[7]*(self.ECONV/self.ACONV)
            atom[8] = atom[8]*(self.ECONV/self.ACONV)
            atom[9] = atom[9]*(self.ECONV/self.ACONV)
            atom[10] = atom[10]*self.ECONV

    def save(self, filename = 'system.json'):

        dictfile = {
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.friction.item(),
                    self.pressure.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'atoms': self.atoms.tolist() 
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

class Symmetry(System):

    def __init__(self):

        self.eta_prm = 0.0
        self.rs_prm = 0.0
        self.rc_prm = 6.0

    def symm2D(self):
        
        for atm in self.atoms:
            sum_ = 0.e0
            for atom in self.atoms:
                if atom[0] != atm[0]:
                    dr = np.sqrt((atm[4]-atom[4])**2+(atm[5]-atom[5])**2+(atm[6]-atom[6])**2)
                    sum_ += np.exp(-self.eta_prm*(dr-self.rs_prm)**2)*self.setSF(dr)
            return sum_

    def setSF(self, dr):

        if dr <= self.rc_prm:
            func = 0.5*(np.cos(np.pi*dr/self.rc_prm)+1.e0)
        else:
            func = 0.e0

        return func

