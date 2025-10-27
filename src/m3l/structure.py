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

class System(Constants):

    description = np.array([])
    cell = np.array([], dtype = np.float64)
    temperature = np.array([], dtype = np.float64)
    temp_friction = np.array(0.0e0, dtype = np.float64)
    press_friction = np.array(0.0e0, dtype = np.float64)
    pressure = np.array([], dtype = np.float64)
    epotential = np.array([], dtype = np.float64)
    ekinetic = np.array([], dtype = np.float64)
    inter_pairs = np.array([], dtype=np.int32)
    sites = np.array([], dtype = np.int32)
    nsites = np.array([], dtype = np.int32)
    volume = np.array([], dtype = np.float64)
    virial = np.array([], dtype=np.float64)
    natom = np.array([], dtype=np.int32)
    nfree = np.array([], dtype=np.int32)
    ma = np.array([], dtype=np.int32)
    atype = np.array([], dtype=np.int32)
    mass = np.array([], dtype=np.float64)
    charge = np.array([], dtype=np.float64)
    ra = np.array([], dtype=np.float64)
    va = np.array([], dtype=np.float64)
    fa = np.array([], dtype=np.float64)
    ea = np.array([], dtype=np.float64)

    def loadSystem(self, filename):

        with open(filename, 'r') as file:
            json_file = json.load(file)
            self.description = json_file['description']
            self.cell = np.array(json_file['cell'], dtype=np.float64)

            atom_list = json_file['atom'] 

        ma = []
        atype = []
        mass =[] 
        charge = []
        ra = []
        va = []
        fa = []
        ea = []
        for atom in atom_list:
            ma.append(atom[0])
            atype.append(atom[1])
            mass.append(atom[2])
            charge.append(atom[3])
            ra.append([atom[4], atom[5], atom[6]])
            va.append([atom[7], atom[8], atom[9]])
            fa.append([atom[10], atom[11], atom[12]])
            ea.append(atom[13])

        self.ma = np.array(ma, dtype=np.int32)
        self.atype = np.array(atype, dtype=np.int32)
        self.mass = np.array(mass, dtype=np.float64)
        self.charge = np.array(charge, dtype=np.float64)
        self.ra = np.array(ra, dtype=np.float64)
        self.va = np.array(va, dtype=np.float64)
        self.fa = np.array(fa, dtype=np.float64)
        self.ea = np.array(ea, dtype=np.float64)

#        self.setMolecules()

        self.setSites()

        self.setNatom(len(atom_list))

        self.setVolume()

        self.setNFree()

        self.setEkinetic()

        self.setTemperature()

        self.setVirial()

        self.setPressure()

        self.setPotential()

    def setSystem(self, description, temperature, pressure, cell, atoms):

        self.description = description
        self.temperature = np.array(temperature, dtype = np.float64)
        self.pressure = np.array(pressure, dtype = np.float64)
        self.virial = np.array(0.0, dtype = np.float64)
        self.cell = np.array(cell, dtype = np.float64)
        self.epotential = np.array(0.0, dtype = np.float64)
        self.temp_friction = np.array(0.0e0, dtype = np.float64)
        self.press_friction = np.array(0.0e0, dtype = np.float64)

        ra = []
        mass = []
        atype = []
        ma = []
        charge = []
        fa = []
        ea = []
        for atom in atoms:
            ma.append(atom[0])
            atype.append(atom[1])
            mass.append(atom[2])
            charge.append(atom[3])
            ra.append([atom[4], atom[5], atom[6]])
            fa.append([0.0, 0.0, 0.0])
            ea.append(0.0)
            
        self.ra = np.array(ra, dtype = np.float64)
        self.mass = np.array(mass, dtype = np.float64)
        self.charge = np.array(charge, dtype = np.float64)
        self.ea = np.array(ea, dtype = np.float64)
        self.fa = np.array(fa, dtype = np.float64)
        self.atype = np.array(atype, dtype = np.int32)
        self.ma = np.array(ma, dtype = np.int32)

        self.setNatom(len(atoms))

        self.setNFree()

        self.setVelocity()

        self.setEkinetic()

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
        for i in range(self.natom):
            sum_ += self.mass[i]*(self.va[i][0]**2+self.va[i][1]**2+self.va[i][2]**2)

        self.ekinetic = np.array(0.5e0*sum_, dtype=np.float64)

        return self.ekinetic

    def setVelocity(self):

        lista = []
        for i in range(self.natom):
            var = np.sqrt(self.nfree*self.temperature/self.mass[i])
            lista.append([var, var, var])

        self.va = np.array(lista, dtype=np.float64)

    def setTemperature(self):

        self.temperature = np.array(2.0*self.ekinetic/self.nfree, dtype=np.float64)

        return self.temperature

    def setVirial(self):

        sum_ = 0.0e0
        for i in range(self.natom):
            for j in range(i+1, self.natom):
                drx = self.ra[j][0]-self.ra[i][0]
                dry = self.ra[j][1]-self.ra[i][1]
                drz = self.ra[j][2]-self.ra[i][2]
                sum_ += self.fa[i][0]*drx+self.fa[i][1]*dry+self.fa[i][2]*drz

        self.virial = np.array(sum_, dtype=np.float64)

        return self.virial

    def setPressure(self):

        self.pressure = np.array((2.0*self.ekinetic+self.virial)/(3.0*self.volume), dtype = np.float64)

        return self.pressure

    def setPotential(self):

        self.epotential = np.array(sum(self.ea), dtype = np.float64)

        return self.epotential

    def setSites(self):

        list_ = []
        list2_ = []
        nx = 0
        for i in range(len(self.atype)):
            for j in range(i+1, len(self.atype)):
                if self.ma[j] is not self.ma[i]:
                    list_.append([i, j])
                    nx = 1 
                else:
                    nx += 1 
                list2_.append(nx)

        self.sites = np.array(list_, dtype=np.int32)
        self.nsites = np.array(list2_, dtype=np.int32) 

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

        self.mass = np.multiply(self.mass, self.MCONV)
        self.charge = np.multiply(self.charge, self.CHG)
        self.ra = np.multiply(self.ra, self.ACONV)
        self.va = np.multiply(self.va, self.ACONV/self.TIMECONV)
        self.fa = np.multiply(self.fa, self.ECONV/self.ACONV)
        self.ea = np.multiply(self.ea, self.ECONV)

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

        self.mass = np.divide(self.mass, self.MCONV)
        self.charge = np.divide(self.charge, self.CHG)
        self.ra = np.divide(self.ra, self.ACONV)
        self.va = np.divide(self.va, self.ACONV/self.TIMECONV)
        self.fa = np.divide(self.fa, self.ECONV/self.ACONV)
        self.ea = np.divide(self.ea, self.ECONV)

    def save(self, filename = 'system.json'):

        atom_list = []
        for i in range(self.natom):
            atom_list.append([self.ma[i].item(),
                             self.atype[i].item(),
                             self.mass[i].item(), 
                             self.charge[i].item(),
                             self.ra[i][0].item(),
                             self.ra[i][1].item(),
                             self.ra[i][2].item(),
                             self.va[i][0].item(),
                             self.va[i][1].item(),
                             self.va[i][2].item(),
                             self.fa[i][0].item(),
                             self.fa[i][1].item(),
                             self.fa[i][2].item(),
                             self.ea[i].item()]
                             )

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
                'atom': atom_list
                }

        outfile = json.dumps(dictfile, indent = 1)
        with open(filename, 'w') as file:
            file.write(outfile)

