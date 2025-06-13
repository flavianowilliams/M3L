import numpy as np
from m3l.structure import System
from m3l.utils import Constants
from m3l import libens, libtherm, libstruct, libforce

class Ensemble(Constants):

    def __init__(self, timestep, force_field, temp_bath = None, press_bath = None, tstat = None, pstat = None, bfactor = None):

        self.timestep = np.array(timestep*self.TIMECONV, dtype = np.float64)
        self.dtimestep = np.array(timestep*self.TIMECONV, dtype = np.float64)
        self.force_field = np.array(force_field, dtype = np.float64)
        self.friction = 0.e0

        if tstat:
            self.tstat = np.array(tstat*self.TIMECONV, dtype = np.float64)
        else:
            self.tstat = np.array(2.0*self.TIMECONV, dtype = np.float64)

        if pstat:
            self.pstat = np.array(pstat*self.TIMECONV, dtype = np.float64)
        else:
            self.pstat = np.array(2.0*self.TIMECONV, dtype = np.float64)

        if temp_bath:
            self.temp_bath = np.array(temp_bath*self.TEMPCONV, dtype = np.float64)
        else:
            self.temp_bath = np.array([], dtype = np.float64)

        if press_bath:
            self.press_bath = np.array(press_bath*self.PCONV, dtype = np.float64)
        else:
            self.press_bath = np.array([], dtype = np.float64)

        if bfactor:
            self.bfactor = np.array(bfactor, dtype = np.float64)
        else:
            self.bfactor = self.BETAFACTOR
    
        self.bfactor = self.bfactor*self.PCONVINV

    def nve_verlet(self):

        print("não valendo...")

    def nvt_verlet(self):

        libens.temp_friction = self.sys.temp_friction
        libens.press_friction = self.sys.press_friction

#definindo coordenadas canonicas

        libstruct.natom = self.sys.natom
        libstruct.nfree = self.nfree
        libstruct.cell = self.sys.cell
        libstruct.atom = self.sys.atom

#-atualizando coordenadas canônicas

        for atom in libstruct.atom:
            atom[6] = atom[6]+0.5e0*self.timestep*atom[9]/atom[1]
            atom[7] = atom[7]+0.5e0*self.timestep*atom[10]/atom[1]
            atom[8] = atom[8]+0.5e0*self.timestep*atom[11]/atom[1]
            atom[3] = atom[3]+atom[6]*self.timestep
            atom[4] = atom[4]+atom[7]*self.timestep
            atom[5] = atom[5]+atom[8]*self.timestep

#-atualizando coordenadas atômicas
        
        libstruct.ccp()

#-calculando campo de força

        libforce.natom = libstruct.natom

        libforce.atom = libstruct.atom
        libforce.nmolecules = len(self.sys.molecule)
        libforce.sites = self.sys.sites
        libforce.nsites = self.sys.nsites
        libforce.molecules = self.sys.molecule
        libforce.nvdw = len(self.force_field)-1
        libforce.rvdw = self.force_field[0, 0]
        libforce.rcoul = self.force_field[0, 1]
        libforce.cell = libstruct.cell
        libforce.params = self.force_field 
        libforce.setforces()

        libstruct.atom = libforce.atom

#-atualizando coordenadas canônicas

        for atom in libstruct.atom:
            atom[6] = atom[6]+0.5e0*self.timestep*atom[9]/atom[1]
            atom[7] = atom[7]+0.5e0*self.timestep*atom[10]/atom[1]
            atom[8] = atom[8]+0.5e0*self.timestep*atom[11]/atom[1]

#-calculando energia cinética

        libtherm.natom = libstruct.natom
        libtherm.atom = libstruct.atom
        libtherm.setekinetic()

#-atualizando coordenadas canônicas

        sigma = 0.5e0*self.nfree*self.temp_bath
        qui = np.sqrt(1.0e0+self.timestep*(sigma/libtherm.ekinetic-1.0e0)/self.tstat)

        for atom in libstruct.atom:
            atom[6] = atom[6]*qui
            atom[7] = atom[7]*qui
            atom[8] = atom[8]*qui

#-calculando energia cinética e temperatura

        libtherm.atom = libstruct.atom
        libtherm.nfree = libstruct.nfree
        libtherm.setekinetic()
        libtherm.settemperature()

#-calculando pressão 

        libtherm.virial = libforce.virial
        libtherm.setpressure()

#-definindo volume 

        libtherm.volume = self.sys.volume

    def npt_hoover(self):

        libens.temp_friction = self.sys.temp_friction
        libens.press_friction = self.sys.press_friction

#definindo coordenadas canonicas

        libstruct.natom = self.sys.natom
        libstruct.nfree = self.nfree
        libstruct.cell = self.sys.cell
        libstruct.atom = self.sys.atom

    def npt_verlet(self):

        libens.temp_friction = self.sys.temp_friction
        libens.press_friction = self.sys.press_friction

#definindo coordenadas canonicas

        libstruct.natom = self.sys.natom
        libstruct.nfree = self.nfree
        libstruct.cell = self.sys.cell
        libstruct.atom = self.sys.atom

#- calculando parametro eta

        libens.bfactor = self.bfactor
        libens.timestep = self.timestep
        libens.press_bath = self.press_bath
        libens.pressure = self.sys.pressure
        libens.pstat = self.pstat

        libens.seteta()
        eta = libens.eta

#-atualizando coordenadas canônicas

#        libens.natom = self.sys.natom
#        libens.atom = atm
#        libens.npt_berendsen_a()
#        libens.npt_berendsen_b()

        for atom in libstruct.atom:
            atom[6] = atom[6]+0.5e0*self.timestep*atom[9]/atom[1]
            atom[7] = atom[7]+0.5e0*self.timestep*atom[10]/atom[1]
            atom[8] = atom[8]+0.5e0*self.timestep*atom[11]/atom[1]
            atom[3] = atom[3]*eta+atom[6]*self.timestep
            atom[4] = atom[4]*eta+atom[7]*self.timestep
            atom[5] = atom[5]*eta+atom[8]*self.timestep

#-atualizando parâmetros de rede e volume da supercélula

        libstruct.eta = eta
        libstruct.setcell()

        libtherm.cell = libstruct.cell
        libtherm.volume = self.sys.volume
        libtherm.eta = eta
        libtherm.setvolumex()

#-atualizando coordenadas atômicas
        
        libstruct.ccp()

#-calculando campo de força

        libforce.natom = libstruct.natom

        libforce.atom = libstruct.atom
        libforce.nmolecules = len(self.sys.molecule)
        libforce.sites = self.sys.sites
        libforce.nsites = self.sys.nsites
        libforce.molecules = self.sys.molecule
        libforce.nvdw = len(self.force_field)-1
        libforce.rvdw = self.force_field[0, 0]
        libforce.rcoul = self.force_field[0, 1]
        libforce.cell = libstruct.cell
        libforce.params = self.force_field 
        libforce.setforces()

        libstruct.atom = libforce.atom

#-atualizando coordenadas canônicas

#        libens.atom = libstruct.atom
#        libens.npt_berendsen_a()

        for atom in libstruct.atom:
            atom[6] = atom[6]+0.5e0*self.timestep*atom[9]/atom[1]
            atom[7] = atom[7]+0.5e0*self.timestep*atom[10]/atom[1]
            atom[8] = atom[8]+0.5e0*self.timestep*atom[11]/atom[1]

#-calculando energia cinética

        libtherm.natom = libstruct.natom
        libtherm.atom = libstruct.atom
        libtherm.setekinetic()

#-atualizando coordenadas canônicas

#        libens.nfree = libstruct.nfree
#        libens.ekinetic = libtherm.ekinetic
#        libens.tstat = self.tstat
#        libens.temp_bath = self.temp_bath
#        libens.npt_berendsen_c()

        sigma = 0.5e0*self.nfree*self.temp_bath
        qui = np.sqrt(1.0e0+self.timestep*(sigma/libtherm.ekinetic-1.0e0)/self.tstat)

        for atom in libstruct.atom:
            atom[6] = atom[6]*qui
            atom[7] = atom[7]*qui
            atom[8] = atom[8]*qui

#-calculando energia cinética e temperatura

        libtherm.atom = libstruct.atom
        libtherm.nfree = libstruct.nfree
        libtherm.setekinetic()
        libtherm.settemperature()

#-calculando pressão 

        libtherm.virial = libforce.virial
        libtherm.setpressure()

    def hook(self, system):

        self.sys = system
        self.nfree = 3*(system.natom-1)

    def hook_output(self):

        description = self.sys.description

        self.sys = System()

        self.sys.description = description
        self.sys.cell = libstruct.cell
        self.sys.volume = libtherm.volume
        self.sys.atom = libstruct.atom
        self.sys.natom = libstruct.natom
        self.sys.sites = libforce.sites
        self.sys.nsites = libforce.nsites
        self.sys.molecule = libforce.molecules
#        
        self.sys.temperature = libtherm.temperature
        self.sys.temp_friction = libens.temp_friction
        self.sys.pressure = libtherm.pressure
        self.sys.press_friction = libens.press_friction
        self.sys.ekinetic = libtherm.ekinetic
        self.sys.epotential = libforce.energy
        self.sys.volume = libtherm.volume
        self.sys.virial = libforce.virial

        return self.sys

    def __call__(self, system):
        self.hook(system)
        if self.temp_bath.size > 0:
            if self.press_bath.size > 0:
                self.npt_verlet()
            else:
                self.nvt_verlet()
        else:
            self.nve_verlet()
        output = self.hook_output()
        return output 

