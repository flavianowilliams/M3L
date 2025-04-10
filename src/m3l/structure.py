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
        
        self.cell = cell
        self.temperature = temperature
        self.pressure = pressure
        self.epotential = epotential
        self.ekinetic = ekinetic
        self.mat = mat
        self.atp = atp
        self.rx = rx
        self.ry = ry
        self.rz = rz
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.fx = fx
        self.fy = fy
        self.fz = fz
        self.ea = ea

    def convertUnitsInv(self):

        self.epotential = np.divide(self.epotential, self.ECONV)
        self.ekinetic = np.divide(self.ekinetic, self.ECONV)
        self.temperature = np.divide(self.temperature, self.TEMPCONV)
        self.pressure = np.divide(self.pressure, self.PCONV)

        self.cell[0] = np.divide(self.cell[0], self.ACONV)
        self.cell[1] = np.divide(self.cell[1], self.ACONV)
        self.cell[2] = np.divide(self.cell[2], self.ACONV)

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

#    def save(self, filename = 'system.json'):
#
#        dictfile = {
#                'cell': self.cell.tolist(),
#                'thermodynamic': [
#                    self.temperature.item(),
#                    self.pressure.item(),
#                    self.epotential.item(),
#                    self.ekinetic.item()
#                    ],
#                'mat': self.mat.tolist(),
#                'atp': self.atp.tolist(),
#                'rx': self.rx.tolist(),
#                'ry': self.ry.tolist(),
#                'rz': self.rz.tolist(),
#                'vx': self.vx.tolist(),
#                'vy': self.vy.tolist(),
#                'vz': self.vz.tolist(),
#                'fx': self.fx.tolist(),
#                'fy': self.fy.tolist(),
#                'fz': self.fz.tolist(),
#                'ea': self.ea.tolist(),
#                }
#
#        outfile = json.dumps(dictfile, indent = 1)
#        with open(filename, 'w') as file:
#            file.write(outfile)

    def save(self, filename = 'system.json'):

        atoms = np.array([[
                self.mat[0],
                self.atp[0],
                self.rx[0],
                self.ry[0],
                self.rz[0],
                self.vx[0],
                self.vy[0],
                self.vz[0],
                self.fx[0],
                self.fy[0],
                self.fz[0],
                self.ea[0]
            ]], dtype = np.float64)

        for i in range(1, len(self.mat)):

            atoms = np.append(atoms, [[
                self.mat[i],
                self.atp[i],
                self.rx[i],
                self.ry[i],
                self.rz[i],
                self.vx[i],
                self.vy[i],
                self.vz[i],
                self.fx[i],
                self.fy[i],
                self.fz[i],
                self.ea[i]
                ]], axis = 0)

        dictfile = {
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.pressure.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'atoms': atoms.tolist()
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

class System(Constants):

    def loadSystem(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.cell = np.array(json_file['cell'], dtype=np.float64)
            self.temperature = np.array(json_file['thermodynamic'][0], dtype=np.float64)
            self.pressure = np.array(json_file['thermodynamic'][1], dtype=np.float64)
            self.epotential = np.array(json_file['thermodynamic'][2], dtype=np.float64)
            self.ekinetic = np.array(json_file['thermodynamic'][3], dtype=np.float64)

            atoms = np.array(json_file['atoms'], dtype=np.float64)

            self.mat = np.array([])
            self.atp = np.array([])
            self.rx = np.array([])
            self.ry = np.array([])
            self.rz = np.array([])
            self.vx = np.array([])
            self.vy = np.array([])
            self.vz = np.array([])
            self.fx = np.array([])
            self.fy = np.array([])
            self.fz = np.array([])
            self.ea = np.array([])

            for atom in atoms:

                self.mat = np.append(self.mat, atom[0])
                self.atp = np.append(self.atp, atom[1])
                self.rx = np.append(self.rx, atom[2])
                self.ry = np.append(self.ry, atom[3])
                self.rz = np.append(self.rz, atom[4])
                self.vx = np.append(self.vx, atom[5])
                self.vy = np.append(self.vy, atom[6])
                self.vz = np.append(self.vz, atom[7])
                self.fx = np.append(self.fx, atom[8])
                self.fy = np.append(self.fy, atom[9])
                self.fz = np.append(self.fz, atom[10])
                self.ea = np.append(self.ea, atom[11])

        self.setSites()

# convertendo datatypes

        self.atp = np.array(self.atp, dtype = np.int32)

#    def loadSystem2(self, filename):
#
#        with open(filename, 'r') as file:
#            json_file = json.load(file)
#            self.cell = np.array(json_file['cell'], dtype=np.float64)
#            self.temperature = np.array(json_file['thermodynamic'][0], dtype=np.float64)
#            self.pressure = np.array(json_file['thermodynamic'][1], dtype=np.float64)
#            self.epotential = np.array(json_file['thermodynamic'][2], dtype=np.float64)
#            self.ekinetic = np.array(json_file['thermodynamic'][3], dtype=np.float64)
#            self.mat = np.array(json_file['mat'], dtype=np.float64)
#            self.atp = np.array(json_file['atp'], dtype=np.int32)
#            self.rx = np.array(json_file['rx'], dtype=np.float64)
#            self.ry = np.array(json_file['ry'], dtype=np.float64)
#            self.rz = np.array(json_file['rz'], dtype=np.float64)
#            self.vx = np.array(json_file['vx'], dtype=np.float64)
#            self.vy = np.array(json_file['vy'], dtype=np.float64)
#            self.vz = np.array(json_file['vz'], dtype=np.float64)
#            self.fx = np.array(json_file['fx'], dtype=np.float64)
#            self.fy = np.array(json_file['fy'], dtype=np.float64)
#            self.fz = np.array(json_file['fz'], dtype=np.float64)
#            self.ea = np.array(json_file['ea'], dtype=np.float64)

    def setSystem(self, temperature, pressure, cell, atoms):

        self.temperature = np.array(temperature, dtype = np.float64)
        self.pressure = np.array(pressure, dtype = np.float64)
        self.cell = np.array(cell, dtype = np.float64)
        self.epotential = np.array(0.0, dtype = np.float64)

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

        self.mat = np.array(mat, dtype = np.float64)
        self.atp = np.array(atp, dtype = np.int32)
        self.rx = np.array(rx, dtype = np.float64)
        self.ry = np.array(ry, dtype = np.float64)
        self.rz = np.array(rz, dtype = np.float64)
        self.vx = np.repeat(self.ekinetic, natoms)
        self.vy = np.repeat(self.ekinetic, natoms)
        self.vz = np.repeat(self.ekinetic, natoms)
        self.fx = np.zeros(natoms)
        self.fy = np.zeros(natoms)
        self.fz = np.zeros(natoms)
        self.ea = np.zeros(natoms)

    def setEkinetic(self, natoms):

        nfree = 3*(natoms-1)
        
        self.ekinetic = np.array(1.5*nfree*self.KB*(self.temperature*self.TEMPCONV), dtype = np.float64)

    def setSites(self):

        self.site = np.array([], dtype = np.int32)

        for item in self.atp:

            if item not in self.site:

                self.site = np.append(self.site, item)

        self.nsite = np.array([], dtype = np.int32)

        for site in self.site:

            nx = 0 

            for item in self.atp:

                if site == item:

                    nx += 1 

            self.nsite = np.append(self.nsite, nx)

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

        self.epotential = np.divide(self.epotential, self.ECONV)
        self.ekinetic = np.divide(self.ekinetic, self.ECONV)
        self.temperature = np.divide(self.temperature, self.TEMPCONV)
        self.pressure = np.divide(self.pressure, self.PCONV)

        self.cell[0] = np.divide(self.cell[0], self.ACONV)
        self.cell[1] = np.divide(self.cell[1], self.ACONV)
        self.cell[2] = np.divide(self.cell[2], self.ACONV)

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

        atoms = np.array([[
                self.mat[0],
                self.atp[0],
                self.rx[0],
                self.ry[0],
                self.rz[0],
                self.vx[0],
                self.vy[0],
                self.vz[0],
                self.fx[0],
                self.fy[0],
                self.fz[0],
                self.ea[0]
            ]], dtype = np.float64)

        for i in range(1, len(self.mat)):

            atoms = np.append(atoms, [[
                self.mat[i],
                self.atp[i],
                self.rx[i],
                self.ry[i],
                self.rz[i],
                self.vx[i],
                self.vy[i],
                self.vz[i],
                self.fx[i],
                self.fy[i],
                self.fz[i],
                self.ea[i]
                ]], axis = 0)

        dictfile = {
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.pressure.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'atoms': atoms.tolist()
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

#    def save(self, filename = 'system.json'):
#
#        dictfile = {
#                'cell': self.cell.tolist(),
#                'thermodynamic': [
#                    self.temperature.item(),
#                    self.pressure.item(),
#                    self.epotential.item(),
#                    self.ekinetic.item()
#                    ],
#                'mat': self.mat.tolist(),
#                'atp': self.atp.tolist(),
#                'rx': self.rx.tolist(),
#                'ry': self.ry.tolist(),
#                'rz': self.rz.tolist(),
#                'vx': self.vx.tolist(),
#                'vy': self.vy.tolist(),
#                'vz': self.vz.tolist(),
#                'fx': self.fx.tolist(),
#                'fy': self.fy.tolist(),
#                'fz': self.fz.tolist(),
#                'ea': self.ea.tolist(),
#                }
#
#        outfile = json.dumps(dictfile, indent = 1)
#        with open(filename, 'w') as file:
#            file.write(outfile)

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

