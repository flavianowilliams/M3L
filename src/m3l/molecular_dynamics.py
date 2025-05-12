import numpy as np
from m3l.structure import Atom, System
from m3l.utils import Constants
from m3l import libs

class Ensemble(Constants):

    def __init__(self, timestep, force_field, temp_bath = None, press_bath = None, tstat = None, pstat = None, bfactor = None):

        self.timestep = np.array(timestep*self.TIMECONV, dtype = np.float64)
        self.dtimestep = np.array(timestep*self.TIMECONV, dtype = np.float64)
        self.force_field = np.array(force_field, dtype = np.float64)
        self.temp_bath = np.array(temp_bath*self.TEMPCONV, dtype = np.float64)
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

        print("não valendo...")

    def npt_verlet(self):

        libs.nfree = self.nfree
        libs.timestep = self.dtimestep
        libs.bfactor = self.bfactor
        libs.tstat = self.tstat
        libs.pstat = self.pstat
        libs.press_bath = self.press_bath
        libs.temp_bath = self.temp_bath
        libs.params = self.force_field 
        libs.nvdw = len(self.force_field)-1
        libs.rvdw = self.force_field[0, 0]
        libs.rcoul = self.force_field[0, 1]

        libs.atom = self.sys.atom
        libs.natom = self.sys.natom
        libs.cell = self.sys.cell
        libs.volume = self.sys.volume
        libs.sites = self.sys.sites
        libs.nsites = self.sys.nsites
        libs.nmolecules = len(self.sys.molecule)
        libs.molecules = self.sys.molecule

        libs.prepare()
        libs.npt_berendsen()

    def hook(self, system):

        self.sys = system
        self.nfree = 3*(system.natom-1)

    def hook_output(self):

        description = self.sys.description

        self.sys = System()

        self.sys.description = description
        self.sys.cell = libs.cell
        self.sys.volume = libs.volume
        self.sys.atom = libs.atom
        self.sys.natom = libs.natom
        self.sys.sites = libs.sites
        self.sys.nsites = libs.nsites
        self.sys.molecule = libs.molecules
        
        self.sys.temperature = libs.temperature
        self.sys.temp_friction = libs.temp_friction
        self.sys.temp_bath = libs.temp_bath
        self.sys.pressure = libs.pressure
        self.sys.press_bath = libs.press_bath
        self.sys.press_friction = libs.press_friction
        self.sys.ekinetic = libs.ekinetic
        self.sys.epotential = libs.energy
        self.sys.volume = libs.volume

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

