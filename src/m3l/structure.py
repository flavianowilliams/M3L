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

        return mass*self.MCONV

    def setZNumber(self, atom):
        
        zatom = None

        if atom == 'H':
            zatom = 1.0 
        elif atom == 'Ar':
            zatom = 18.0

        return zatom

class System2(Constants):

    def __init__(self, cell, temperature, pressure, epotential, ekinetic, mat, atp, rx, ry, rz, vx, vy, vz, fx, fy, fz, ea):
        
        self.cell = np.array(cell)
        self.temperature = np.array(temperature)
        self.pressure = np.array(pressure)
        self.epotential = np.array(epotential)
        self.ekinetic = np.array(ekinetic)
        self.mat = np.array(mat)
        self.atp = np.array(atp)
        self.rx = np.array(rx)
        self.ry = np.array(ry)
        self.rz = np.array(rz)
        self.vx = np.array(vx)
        self.vy = np.array(vy)
        self.vz = np.array(vz)
        self.fx = np.array(fx)
        self.fy = np.array(fy)
        self.fz = np.array(fz)
        self.ea = np.array(ea)

    def convertUnitsInv(self):

        self.epotential = self.epotential/self.ECONV
        self.ekinetic = self.ekinetic/self.ECONV
        self.temperature = self.temperature/self.TEMPCONV
        self.pressure = self.pressure/self.PCONV

        self.cell[0] = self.cell[0]/self.ACONV
        self.cell[1] = self.cell[1]/self.ACONV
        self.cell[2] = self.cell[2]/self.ACONV

        self.mat = np.divide(self.mat, self.MCONV)
        self.rx = np.divide(self.rx, self.ACONV)
        self.ry = np.divide(self.ry, self.ACONV)
        self.rz = np.divide(self.rz, self.ACONV)
        self.vx = np.divide(self.vx, self.ACONV/self.TIMECONV)
        self.vy = np.divide(self.vy, self.ACONV/self.TIMECONV)
        self.vz = np.divide(self.vz, self.ACONV/self.TIMECONV)
        self.fx = np.divide(self.fx, self.ECONV/self.ACONV)
        self.fy = np.divide(self.fy, self.ECONV/self.ACONV)
        self.fz = np.divide(self.fz, self.ECONV/self.ACONV)
        self.ea = np.divide(self.ea, self.ECONV)

    def save(self, filename = 'system.json'):

        dictfile = {
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.pressure.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'mat': self.mat.tolist(),
                'atp': self.atp.tolist(),
                'rx': self.rx.tolist(),
                'ry': self.ry.tolist(),
                'rz': self.rz.tolist(),
                'vx': self.vx.tolist(),
                'vy': self.vy.tolist(),
                'vz': self.vz.tolist(),
                'fx': self.fx.tolist(),
                'fy': self.fy.tolist(),
                'fz': self.fz.tolist(),
                'ea': self.ea.tolist(),
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

class System(Constants):

    def loadSystem(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.cell = np.array(json_file['cell'], dtype=np.float32)
            self.temperature = np.array(json_file['thermodynamic'][0], dtype=np.float32)
            self.pressure = np.array(json_file['thermodynamic'][1], dtype=np.float32)
            self.epotential = np.array(json_file['thermodynamic'][2], dtype=np.float32)
            self.ekinetic = np.array(json_file['thermodynamic'][3], dtype=np.float32)
            self.mat = np.array(json_file['mat'], dtype=np.float32)
            self.atp = np.array(json_file['atp'], dtype=np.int32)
            self.rx = np.array(json_file['rx'], dtype=np.float32)
            self.ry = np.array(json_file['ry'], dtype=np.float32)
            self.rz = np.array(json_file['rz'], dtype=np.float32)
            self.vx = np.array(json_file['vx'], dtype=np.float32)
            self.vy = np.array(json_file['vy'], dtype=np.float32)
            self.vz = np.array(json_file['vz'], dtype=np.float32)
            self.fx = np.array(json_file['fx'], dtype=np.float32)
            self.fy = np.array(json_file['fy'], dtype=np.float32)
            self.fz = np.array(json_file['fz'], dtype=np.float32)
            self.ea = np.array(json_file['ea'], dtype=np.float32)

    def setSystem(self, temperature, pressure, cell, atoms):

        self.temperature = np.array(temperature)
        self.pressure = np.array(pressure)
        self.cell = np.array(cell)
        self.epotential = np.array(0.0)

        mat = []
        atp = []
        rx = []
        ry = []
        rz = []
        for atom in atoms:
           mat.append(atom[0])
           atp.append(atom[1])
           rx.append(atom[2])
           ry.append(atom[3])
           rz.append(atom[4])

        natoms = len(atoms)

        self.setEkinetic(natoms)

        self.mat = np.array(mat, dtype = np.float32)
        self.atp = np.array(atp, dtype = np.int32)
        self.rx = np.array(rx, dtype = np.float32)
        self.ry = np.array(ry, dtype = np.float32)
        self.rz = np.array(rz, dtype = np.float32)
        self.vx = np.repeat(self.ekinetic, natoms)
        self.vy = np.repeat(self.ekinetic, natoms)
        self.vz = np.repeat(self.ekinetic, natoms)
        self.fx = np.zeros(natoms)
        self.fy = np.zeros(natoms)
        self.fz = np.zeros(natoms)
        self.ea = np.zeros(natoms)

    def setEkinetic(self, natoms):

        nfree = 3*(natoms-1)
        
        self.ekinetic = np.array(1.5*nfree*self.KB*(self.temperature*self.TEMPCONV), dtype = np.float32)

    def convertUnits(self):

        self.epotential = self.epotential*self.ECONV
        self.ekinetic = self.ekinetic*self.ECONV
        self.temperature = self.temperature*self.TEMPCONV

        self.cell[0] = self.cell[0]*self.ACONV
        self.cell[1] = self.cell[1]*self.ACONV
        self.cell[2] = self.cell[2]*self.ACONV

        self.mat = np.multiply(self.mat, self.MCONV)
        self.rx = np.multiply(self.rx, self.ACONV)
        self.ry = np.multiply(self.ry, self.ACONV)
        self.rz = np.multiply(self.rz, self.ACONV)
        self.vx = np.multiply(self.vx, self.ACONV/self.TIMECONV)
        self.vy = np.multiply(self.vy, self.ACONV/self.TIMECONV)
        self.vz = np.multiply(self.vz, self.ACONV/self.TIMECONV)
        self.fx = np.multiply(self.fx, self.ECONV/self.ACONV)
        self.fy = np.multiply(self.fy, self.ECONV/self.ACONV)
        self.fz = np.multiply(self.fz, self.ECONV/self.ACONV)
        self.ea = np.multiply(self.ea, self.ECONV)

    def convertUnitsInv(self):

        self.epotential = self.epotential/self.ECONV
        self.ekinetic = self.ekinetic/self.ECONV
        self.temperature = self.temperature/self.TEMPCONV
        self.pressure = self.pressure/self.PCONV

        self.cell[0] = self.cell[0]/self.ACONV
        self.cell[1] = self.cell[1]/self.ACONV
        self.cell[2] = self.cell[2]/self.ACONV

        self.mat = np.divide(self.mat, self.MCONV)
        self.rx = np.divide(self.rx, self.ACONV)
        self.ry = np.divide(self.ry, self.ACONV)
        self.rz = np.divide(self.rz, self.ACONV)
        self.vx = np.divide(self.vx, self.ACONV/self.TIMECONV)
        self.vy = np.divide(self.vy, self.ACONV/self.TIMECONV)
        self.vz = np.divide(self.vz, self.ACONV/self.TIMECONV)
        self.fx = np.divide(self.fx, self.ECONV/self.ACONV)
        self.fy = np.divide(self.fy, self.ECONV/self.ACONV)
        self.fz = np.divide(self.fz, self.ECONV/self.ACONV)
        self.ea = np.divide(self.ea, self.ECONV)

    def save(self, filename = 'system.json'):

        dictfile = {
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.pressure.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'mat': self.mat.tolist(),
                'atp': self.atp.tolist(),
                'rx': self.rx.tolist(),
                'ry': self.ry.tolist(),
                'rz': self.rz.tolist(),
                'vx': self.vx.tolist(),
                'vy': self.vy.tolist(),
                'vz': self.vz.tolist(),
                'fx': self.fx.tolist(),
                'fy': self.fy.tolist(),
                'fz': self.fz.tolist(),
                'ea': self.ea.tolist(),
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

