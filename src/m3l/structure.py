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

#class System2(Constants):
#
#    def __init__(self, description, cell, volume, temperature, temp_friction, pressure, press_friction, epotential, ekinetic, nsites, sites, molecules, mat, atp, chg, rx, ry, rz, vx, vy, vz, fx, fy, fz, ea):
#        
#        self.description = description
#        self.cell = cell
#        self.volume = volume
#        self.temperature = temperature
#        self.temp_friction = temp_friction
#        self.pressure = pressure
#        self.press_friction = press_friction
#        self.epotential = epotential
#        self.ekinetic = ekinetic
#        self.nsites = nsites
#        self.sites = sites
#        self.molecules = molecules
#        self.mat = mat
#        self.atp = atp
#        self.chg = chg
#        self.rx = rx
#        self.ry = ry
#        self.rz = rz
#        self.vx = vx
#        self.vy = vy
#        self.vz = vz
#        self.fx = fx
#        self.fy = fy
#        self.fz = fz
#        self.ea = ea
#
#    def convertUnitsInv(self):
#
#        self.epotential = np.divide(self.epotential, self.ECONV)
#        self.ekinetic = np.divide(self.ekinetic, self.ECONV)
#        self.temperature = np.divide(self.temperature, self.TEMPCONV)
#        self.pressure = np.divide(self.pressure, self.PCONV)
#        self.temp_friction = np.divide(self.temp_friction, self.TEMPCONV/self.TIMECONV)
#        self.press_friction = np.divide(self.press_friction, self.PCONV/self.TIMECONV)
#
#
#        self.cell[0] = np.divide(self.cell[0], self.ACONV)
#        self.cell[1] = np.divide(self.cell[1], self.ACONV)
#        self.cell[2] = np.divide(self.cell[2], self.ACONV)
#        self.volume = np.divide(self.volume, self.ACONV**3)
#
#        self.mat = np.divide(self.mat, self.MCONV)
#        self.rx = np.divide(self.rx, self.ACONV)
#        self.ry = np.divide(self.ry, self.ACONV)
#        self.rz = np.divide(self.rz, self.ACONV)
#        self.vx = np.divide(self.vx, self.ACONV/self.TIMECONV)
#        self.vy = np.divide(self.vy, self.ACONV/self.TIMECONV)
#        self.vz = np.divide(self.vz, self.ACONV/self.TIMECONV)
#        self.fx = np.divide(self.fx, self.ECONV/self.ACONV)
#        self.fy = np.divide(self.fy, self.ECONV/self.ACONV)
#        self.fz = np.divide(self.fz, self.ECONV/self.ACONV)
#        self.ea = np.divide(self.ea, self.ECONV)
#
#    def save(self, filename = 'system.json'):
#
#        atoms = np.array([[
#                self.atp[0],
#                self.mat[0],
#                self.chg[0],
#                self.rx[0],
#                self.ry[0],
#                self.rz[0],
#                self.vx[0],
#                self.vy[0],
#                self.vz[0],
#                self.fx[0],
#                self.fy[0],
#                self.fz[0],
#                self.ea[0]
#            ]], dtype = np.float64)
#
#        for i in range(1, len(self.mat)):
#
#            atoms = np.append(atoms, [[
#                self.atp[i],
#                self.mat[i],
#                self.chg[i],
#                self.rx[i],
#                self.ry[i],
#                self.rz[i],
#                self.vx[i],
#                self.vy[i],
#                self.vz[i],
#                self.fx[i],
#                self.fy[i],
#                self.fz[i],
#                self.ea[i]
#                ]], axis = 0)
#
#        dictfile = {
#                'description': self.description,
#                'cell': self.cell.tolist(),
#                'molecules': self.molecules.tolist(),
#                'thermodynamic': [
#                    self.temperature.item(),
#                    self.temp_friction.item(),
#                    self.pressure.item(),
#                    self.press_friction.item(),
#                    self.epotential.item(),
#                    self.ekinetic.item()
#                    ],
#                'atoms': atoms.tolist()
#                }
#
#        outfile = json.dumps(dictfile, indent = 1)
#        with open(filename, 'w') as file:
#            file.write(outfile)

class System(Constants):

    description = np.array([])
    cell = np.array([], dtype = np.float64)
    temperature = np.array([], dtype = np.float64)
    temp_friction = np.array([], dtype = np.float64)
    pressure = np.array([], dtype = np.float64)
    press_friction = np.array([], dtype = np.float64)
    epotential = np.array([], dtype = np.float64)
    ekinetic = np.array([], dtype = np.float64)
    molecule = np.array([], dtype = np.int32)
    sites = np.array([], dtype = np.int32)
    nsites = np.array([], dtype = np.int32)
    atom = np.array([], dtype = np.float64)
    volume = np.array([], dtype = np.float64)
    atp = np.array([], dtype = np.float64)
    mat = np.array([], dtype = np.float64)
    chg = np.array([], dtype = np.float64)
    rx = np.array([], dtype = np.float64)
    ry = np.array([], dtype = np.float64)
    rz = np.array([], dtype = np.float64)
    vx = np.array([], dtype = np.float64)
    vy = np.array([], dtype = np.float64)
    vz = np.array([], dtype = np.float64)
    fx = np.array([], dtype = np.float64)
    fy = np.array([], dtype = np.float64)
    fz = np.array([], dtype = np.float64)
    ea = np.array([], dtype = np.float64)

    def loadSystem(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.description = json_file['description']
            self.cell = np.array(json_file['cell'], dtype=np.float64)
            self.temperature = np.array(json_file['thermodynamic'][0], dtype=np.float64)
            self.temp_friction = np.array(json_file['thermodynamic'][1], dtype=np.float64)
            self.pressure = np.array(json_file['thermodynamic'][2], dtype=np.float64)
            self.press_friction = np.array(json_file['thermodynamic'][3], dtype=np.float64)
            self.epotential = np.array(json_file['thermodynamic'][4], dtype=np.float64)
            self.ekinetic = np.array(json_file['thermodynamic'][5], dtype=np.float64)
            self.molecule = np.array(json_file['molecule'], dtype=np.int32)

            self.atom = np.array(json_file['atom'], dtype=np.float64)

            self.atp = np.array([])
            self.mat = np.array([])
            self.chg = np.array([])
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

            for atom in self.atom:

                self.atp = np.append(self.atp, atom[0])
                self.mat = np.append(self.mat, atom[1])
                self.chg = np.append(self.chg, atom[2])
                self.rx = np.append(self.rx, atom[3])
                self.ry = np.append(self.ry, atom[4])
                self.rz = np.append(self.rz, atom[5])
                self.vx = np.append(self.vx, atom[6])
                self.vy = np.append(self.vy, atom[7])
                self.vz = np.append(self.vz, atom[8])
                self.fx = np.append(self.fx, atom[9])
                self.fy = np.append(self.fy, atom[10])
                self.fz = np.append(self.fz, atom[11])
                self.ea = np.append(self.ea, atom[12])

        self.setSites()

        self.setNatom(len(self.atom))

        self.setVolume()

# convertendo datatypes

        self.atp = np.array(self.atp, dtype = np.int32)

    def setSystem(self, description, temperature, pressure, cell, molecules, atoms):

        self.description = description
        self.temperature = np.array(temperature, dtype = np.float64)
        self.temp_friction = np.array(0.0, dtype = np.float64)
        self.pressure = np.array(pressure, dtype = np.float64)
        self.press_friction = np.array(0.0, dtype = np.float64)
        self.cell = np.array(cell, dtype = np.float64)
        self.epotential = np.array(0.0, dtype = np.float64)
        self.molecule = np.array(molecule, np.int32)

        self.atom = np.array([[]], dtype = np.float64)

        atp = []
        mat = []
        chg = []
        rx = []
        ry = []
        rz = []
        for atom in atoms:
           atp.append(atom[0])
           chg.append(atom[1])
           mat.append(atom[2])
           rx.append(atom[3])
           ry.append(atom[4])
           rz.append(atom[5])

        natoms = len(atoms)

        self.setNatom(natoms)

        self.setEkinetic(natoms)

        self.atp = np.array(atp, dtype = np.int32)
        self.mat = np.array(mat, dtype = np.float64)
        self.chg = np.array(chg, dtype = np.float64)
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

    def setNatom(self, natoms):

        self.natom = natoms

        return self.natom

    def setVolume(self):

        self.volume = np.prod(self.cell)

    def setEkinetic(self, natoms):

        nfree = 3*(natoms-1)
        
        self.ekinetic = np.array(1.5*nfree*self.KB*(self.temperature*self.TEMPCONV), dtype = np.float64)

    def setSites(self):

#        self.sites = np.array([], dtype = np.int32)

        list_ = []

        for item in self.atp:

            if item not in list_:

                list_.append(item)

        for site in list_:

            nx = 0 

            for item in self.atp:

                if site == item:

                    nx += 1 

            if self.sites.size == 0:

                self.sites = np.array([[site, nx]], dtype = np.int32)

            else:

                self.sites = np.append([self.sites], [site, nx], axis = 0)

        self.nsites = len(self.sites)

    def convertUnits(self):

        self.epotential = self.epotential*self.ECONV
        self.ekinetic = self.ekinetic*self.ECONV
        self.temperature = self.temperature*self.TEMPCONV
        self.temp_friction = np.multiply(self.temp_friction, self.TEMPCONV/self.TIMECONV)
        self.press_friction = np.multiply(self.press_friction, self.PCONV/self.TIMECONV)

        self.cell[0] = self.cell[0]*self.ACONV
        self.cell[1] = self.cell[1]*self.ACONV
        self.cell[2] = self.cell[2]*self.ACONV
        self.volume = np.multiply(self.volume, self.ACONV**3)

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

        for atom in self.atom:
            atom[1] = np.multiply(atom[1], self.MCONV)
            atom[3] = np.multiply(atom[3], self.ACONV)
            atom[4] = np.multiply(atom[4], self.ACONV)
            atom[5] = np.multiply(atom[5], self.ACONV)
            atom[6] = np.multiply(atom[6], self.ACONV/self.TIMECONV)
            atom[7] = np.multiply(atom[7], self.ACONV/self.TIMECONV)
            atom[8] = np.multiply(atom[8], self.ACONV/self.TIMECONV)
            atom[9] = np.multiply(atom[9], self.ECONV/self.TIMECONV)
            atom[10] = np.multiply(atom[10], self.ECONV/self.TIMECONV)
            atom[11] = np.multiply(atom[11], self.ECONV/self.TIMECONV)
            atom[12] = np.multiply(atom[12], self.ECONV)

    def convertUnitsInv(self):

        self.epotential = np.divide(self.epotential, self.ECONV)
        self.ekinetic = np.divide(self.ekinetic, self.ECONV)
        self.temperature = np.divide(self.temperature, self.TEMPCONV)
        self.pressure = np.divide(self.pressure, self.PCONV)
        self.temp_friction = np.divide(self.temp_friction, self.TEMPCONV/self.TIMECONV)
        self.press_friction = np.divide(self.press_friction, self.PCONV/self.TIMECONV)


        self.cell[0] = np.divide(self.cell[0], self.ACONV)
        self.cell[1] = np.divide(self.cell[1], self.ACONV)
        self.cell[2] = np.divide(self.cell[2], self.ACONV)
        self.volume = np.divide(self.volume, self.ACONV**3)

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

        for atom in self.atom:
            atom[1] = np.divide(atom[1], self.MCONV)
            atom[3] = np.divide(atom[3], self.ACONV)
            atom[4] = np.divide(atom[4], self.ACONV)
            atom[5] = np.divide(atom[5], self.ACONV)
            atom[6] = np.divide(atom[6], self.ACONV/self.TIMECONV)
            atom[7] = np.divide(atom[7], self.ACONV/self.TIMECONV)
            atom[8] = np.divide(atom[8], self.ACONV/self.TIMECONV)
            atom[9] = np.divide(atom[9], self.ECONV/self.TIMECONV)
            atom[10] = np.divide(atom[10], self.ECONV/self.TIMECONV)
            atom[11] = np.divide(atom[11], self.ECONV/self.TIMECONV)
            atom[12] = np.divide(atom[12], self.ECONV)

    def save(self, filename = 'system.json'):

        atoms = np.array([[
                self.atp[0],
                self.mat[0],
                self.chg[0],
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
                self.atp[i],
                self.mat[i],
                self.chg[i],
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
                'description': self.description,
                'cell': self.cell.tolist(),
                'thermodynamic': [
                    self.temperature.item(),
                    self.temp_friction.item(),
                    self.pressure.item(),
                    self.press_friction.item(),
                    self.epotential.item(),
                    self.ekinetic.item()
                    ],
                'molecule': self.molecule.tolist(),
                'atom': atoms.tolist()
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

