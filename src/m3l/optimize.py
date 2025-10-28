from m3l.structure import System
from m3l.structure import Constants
from m3l import libstruct, libforce, libens, libtherm
import numpy as np

class Optimize(Constants):

    def __init__(self, force_field, learning_rate = 1.0e-3):

        self.force_field = force_field
        self.learning_rate = learning_rate

        self.setForceField(force_field)

    def setForceField(self, force_field):

        self.nvdw = np.array(force_field[0], dtype=np.int32)
        self.rvdw = np.array(force_field[1], dtype=np.float64)
        self.rcoul = np.array(force_field[2], dtype=np.float64)
        self.ijinter = np.array(force_field[3][0], dtype=np.float64)
        self.prmvdw = np.array(force_field[3][1], dtype=np.float64)

    def SD(self):

#-calculando campo de força

        self.setForces()

        for i in range(self.sys.natom):
            self.sys.ra[i][0] = self.sys.ra[i][0]+self.learning_rate*self.sys.fa[i][0]
            self.sys.ra[i][1] = self.sys.ra[i][1]+self.learning_rate*self.sys.fa[i][1]
            self.sys.ra[i][2] = self.sys.ra[i][2]+self.learning_rate*self.sys.fa[i][2]

#-calculando velocidades

        self.setVelocity()

#-calculando energia cinética

        libtherm.va = self.sys.va 
        libtherm.setekinetic()

#-calculando pressão 

        libtherm.virial = libforce.virial
        libtherm.cell = self.sys.cell
        libtherm.setvolume()
        libtherm.setpressure()
    
#-aplicando condicoes de contorno periodicas
        
        libstruct.natom = self.sys.natom
        libtherm.nfree = self.sys.nfree 
        libstruct.cell = self.sys.cell
        libstruct.ra = self.sys.ra
        libstruct.ccp()

#-atribuindo variaveis canonicas

        libtherm.setvolume()
        libens.mass = self.sys.mass
        libens.ra = self.sys.ra 
        libens.va = self.sys.va 

    def setForces(self):

        libforce.atype = self.sys.atype
        libforce.charge = self.sys.charge
        libforce.ra = self.sys.ra 
        libforce.fa = self.sys.fa 
        libforce.ea = self.sys.ea 
        libforce.sites = self.sys.sites
        libforce.nsites = self.sys.nsites
        libforce.ma = self.sys.ma
        libforce.cell = self.sys.cell
        libforce.nvdw = self.nvdw 
        libforce.rvdw = self.rvdw
        libforce.rcoul = self.rcoul
        libforce.ijinter = self.ijinter
        libforce.prmvdw = self.prmvdw
        libforce.setforces()

    def setVelocity(self):

        lista = []
        for i in range(len(self.sys.va)):
            var = np.sqrt(self.sys.nfree*self.sys.temperature/self.sys.mass[i])
            lista.append([var, var, var])

        return np.array(lista, dtype=np.float64)

    def hook(self, system):

        self.sys = system
        
        self.SD()

    def hook_output(self):

        description = "Optimized system by gradient descent algorithm"
 
        self.sys = System()

        self.sys.description = description
        self.sys.cell = libstruct.cell
        self.sys.volume = libtherm.volume
        self.sys.natom = libstruct.natom
        self.sys.nfree = libtherm.nfree
        self.sys.sites = libforce.sites
        self.sys.nsites = libforce.nsites
        self.sys.atype = libforce.atype
        self.sys.charge = libforce.charge
        self.sys.mass = libens.mass
        self.sys.ra = libens.ra
        self.sys.va = libens.va
        self.sys.ma = libforce.ma
        self.sys.fa = libforce.fa 
        self.sys.ea = libforce.ea 
        
        self.sys.temperature = libtherm.temperature
        self.sys.press_friction = libens.press_friction
        self.sys.temp_friction = libens.temp_friction
        self.sys.pressure = libtherm.pressure
        self.sys.ekinetic = libtherm.ekinetic
        self.sys.epotential = libforce.energy
        self.sys.virial = libforce.virial
 
        return self.sys

    def __call__(self, system):
        self.hook(system)
        output = self.hook_output()
        return output 

