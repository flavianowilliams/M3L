import numpy as np
import json

from numpy._core.numeric import dtype
from numpy._core.numerictypes import int32
from m3l.utils import Constants

class Atom(Constants):

    def setMass(self, zatom):

        mass = 0.0

        if int(zatom) == 1:
            mass = 1.008
        elif int(zatom) == 18:
            mass = 39.948

        return mass*self.MCONV

    def setZNumber(self, atom):
        
        zatom = None

        if atom == 'H':
            zatom = 1.0 
        elif atom == 'Ar':
            zatom = 18.0

        return zatom

class System(Constants):

    description = np.array([])
    cell = np.array([], dtype = np.float64)
    temperature = np.array([], dtype = np.float64)
    temp_friction = np.array(0.0e0, dtype = np.float64)
    press_friction = np.array(0.0e0, dtype = np.float64)
    pressure = np.array([], dtype = np.float64)
    epotential = np.array([], dtype = np.float64)
    ekinetic = np.array([], dtype = np.float64)
    molecule = np.array([], dtype = np.int32)
    sites = np.array([], dtype = np.int32)
    nsites = np.array([], dtype = np.int32)
    atom = np.array([], dtype = np.float64)
    volume = np.array([], dtype = np.float64)
    virial = np.array([], dtype=np.float64)
    natom = np.array([], dtype=np.int32)
    nfree = np.array([], dtype=np.int32)

    def loadSystem(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.description = json_file['description']
            self.cell = np.array(json_file['cell'], dtype=np.float64)
            self.molecule = np.array(json_file['molecule'], dtype=np.int32)

            self.atom = np.array(json_file['atom'], dtype=np.float64)

        self.setSites()

        self.setNatom(len(self.atom))

        self.setVolume()

        self.setNFree()

        self.setEkinetic()

        self.setTemperature()

        self.setVirial()

        self.setPressure()

        self.setPotential()

    def setSystem(self, description, temperature, pressure, cell, molecule, atoms):

        self.description = description
        self.temperature = np.array(temperature, dtype = np.float64)
        self.pressure = np.array(pressure, dtype = np.float64)
        self.virial = np.array(0.0, dtype = np.float64)
        self.cell = np.array(cell, dtype = np.float64)
        self.epotential = np.array(0.0, dtype = np.float64)
        self.molecule = np.array(molecule, np.int32)

        self.atom = np.array([[]], dtype = np.float64)

    def setNatom(self, natoms):

        self.natom = np.array(natoms, dtype=np.int32)

        return self.natom

    def setNFree(self):

        self.nfree = np.array(3*(self.natom-1), dtype=np.int32)

        return self.nfree

    def setVolume(self):

        self.volume = np.prod(self.cell)

        return self.volume

    def setEkinetic(self):

        sum_ = 0.0e0
        for atom in self.atom:
            sum_ += atom[1]*(atom[6]**2+atom[7]**2+atom[8]**2)

        self.ekinetic = np.array(0.5e0*sum_, dtype=np.float64)

        return self.ekinetic

    def setTemperature(self):

        self.temperature = np.array(2.0*self.ekinetic/self.nfree, dtype=np.float64)

        return self.temperature

    def setVirial(self):

        sum_ = 0.0e0
        for i in range(len(self.atom)):
            for j in range(i+1, len(self.atom)):
                drx = self.atom[j][3]-self.atom[i][3]
                dry = self.atom[j][4]-self.atom[i][4]
                drz = self.atom[j][5]-self.atom[i][5]
                sum_ += self.atom[i][9]*drx+self.atom[i][10]*dry+self.atom[i][11]*drz

        self.virial = np.array(sum_, dtype=np.float64)

        return self.virial

    def setPressure(self):

        self.pressure = np.array((2.0*self.ekinetic+self.virial)/(3.0*self.volume), dtype = np.float64)

        return self.pressure

    def setPotential(self):

        sum_ = 0.0e0
        for atom in self.atom:
            sum_ += atom[12]

        self.epotential = np.array(sum_, dtype = np.float64)

        return self.epotential

    def setSites(self):

        list_ = []

        for item in self.atom:

            if item[0] not in list_:

                list_.append(item[0])

        for site in list_:

            nx = 0 

            for item in self.atom:

                if site == item[0]:

                    nx += 1 

            if self.sites.size == 0:

                self.sites = np.array([[site, nx]], dtype = np.int32)

            else:

                self.sites = np.append([self.sites], [site, nx], axis = 0)

        self.nsites = len(self.sites)

    def convertUnits(self):

        self.epotential = np.multiply(self.epotential, self.ECONV)
        self.ekinetic = np.multiply(self.ekinetic, self.ECONV) 
        self.temperature = np.multiply(self.temperature, self.TEMPCONV) 
        self.pressure = np.multiply(self.pressure, self.PCONV)
        self.virial = np.multiply(self.virial, self.ECONV) 
        self.temp_friction = np.multiply(self.temp_friction, self.TEMPCONV/self.TIMECONV)
        self.press_friction = np.multiply(self.press_friction, self.PCONV/self.TIMECONV)

        self.cell[0] = self.cell[0]*self.ACONV
        self.cell[1] = self.cell[1]*self.ACONV
        self.cell[2] = self.cell[2]*self.ACONV
        self.volume = np.multiply(self.volume, self.ACONV**3)

        for atom in self.atom:
            atom[1] = np.multiply(atom[1], self.MCONV)
            atom[3] = np.multiply(atom[3], self.ACONV)
            atom[4] = np.multiply(atom[4], self.ACONV)
            atom[5] = np.multiply(atom[5], self.ACONV)
            atom[6] = np.multiply(atom[6], self.ACONV/self.TIMECONV)
            atom[7] = np.multiply(atom[7], self.ACONV/self.TIMECONV)
            atom[8] = np.multiply(atom[8], self.ACONV/self.TIMECONV)
            atom[9] = np.multiply(atom[9], self.ECONV/self.ACONV)
            atom[10] = np.multiply(atom[10], self.ECONV/self.ACONV)
            atom[11] = np.multiply(atom[11], self.ECONV/self.ACONV)
            atom[12] = np.multiply(atom[12], self.ECONV)

    def convertUnitsInv(self):

        self.epotential = np.divide(self.epotential, self.ECONV)
        self.ekinetic = np.divide(self.ekinetic, self.ECONV)
        self.temperature = np.divide(self.temperature, self.TEMPCONV)
        self.pressure = np.divide(self.pressure, self.PCONV)
        self.virial = np.divide(self.virial, self.ECONV) 
        self.temp_friction = np.divide(self.temp_friction, self.TEMPCONV/self.TIMECONV)
        self.press_friction = np.divide(self.press_friction, self.PCONV/self.TIMECONV)


        self.cell[0] = np.divide(self.cell[0], self.ACONV)
        self.cell[1] = np.divide(self.cell[1], self.ACONV)
        self.cell[2] = np.divide(self.cell[2], self.ACONV)
        self.volume = np.divide(self.volume, self.ACONV**3)

        for atom in self.atom:
            atom[1] = np.divide(atom[1], self.MCONV)
            atom[3] = np.divide(atom[3], self.ACONV)
            atom[4] = np.divide(atom[4], self.ACONV)
            atom[5] = np.divide(atom[5], self.ACONV)
            atom[6] = np.divide(atom[6], self.ACONV/self.TIMECONV)
            atom[7] = np.divide(atom[7], self.ACONV/self.TIMECONV)
            atom[8] = np.divide(atom[8], self.ACONV/self.TIMECONV)
            atom[9] = np.divide(atom[9], self.ECONV/self.ACONV)
            atom[10] = np.divide(atom[10], self.ECONV/self.ACONV)
            atom[11] = np.divide(atom[11], self.ECONV/self.ACONV)
            atom[12] = np.divide(atom[12], self.ECONV)

    def save(self, filename = 'system.json'):

        dictfile = {
                'description': self.description,
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.temp_friction.item(),
                    self.pressure.item(),
                    self.virial.item(),
                    self.press_friction.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'molecule': self.molecule.tolist(),
                'atom': self.atom.tolist()
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

#class Symmetry(System):
#
#    def __init__(self):
#
#        self.eta_prm = 0.0
#        self.rs_prm = 0.0
#        self.rc_prm = 6.0
#
#    def symm2D(self):
#        
#        for atm in self.atoms:
#            sum_ = 0.e0
#            for atom in self.atoms:
#                if atom[0] != atm[0]:
#                    dr = np.sqrt((atm[4]-atom[4])**2+(atm[5]-atom[5])**2+(atm[6]-atom[6])**2)
#                    sum_ += np.exp(-self.eta_prm*(dr-self.rs_prm)**2)*self.setSF(dr)
#            return sum_
#
#    def setSF(self, dr):
#
#        if dr <= self.rc_prm:
#            func = 0.5*(np.cos(np.pi*dr/self.rc_prm)+1.e0)
#        else:
#            func = 0.e0
#
#        return func

