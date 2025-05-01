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
        #        libs.nfree = self.nfree
        #        libs.temp_bath = self.temp_bath
        #        libs.natom = len(self.mat)
        #        libs.timestep = self.dtimestep
        #        libs.tstat = self.tstat
        #        libs.cell = self.cell
        #        libs.sites = self.sites
        #        libs.nsites = self.nsites
        #        libs.nmolecules = len(self.molecules)
        #        libs.molecules = self.molecules
        #
        #        libs.params = self.force_field 
        #        libs.nvdw = len(self.force_field)-1
        #        libs.rvdw = self.force_field[0, 0]
        #
        #        libs.mass = self.mat
        #        libs.rx = self.rx
        #        libs.ry = self.ry
        #        libs.rz = self.rz
        #        libs.vx = self.vx
        #        libs.vy = self.vy
        #        libs.vz = self.vz
        #        libs.fx = self.fx
        #        libs.fy = self.fy
        #        libs.fz = self.fz
        #        libs.ea = self.ea
        #        libs.atp = self.atp
        #
        #        libs.prepare()
        #        libs.nvt_berendsen()
        #
        #        self.cell = libs.cell
        #        
        #        self.rx = libs.rx
        #        self.ry = libs.ry
        #        self.rz = libs.rz
        #        self.vx = libs.vx
        #        self.vy = libs.vy
        #        self.vz = libs.vz
        #        self.fx = libs.fx
        #        self.fy = libs.fy
        #        self.fz = libs.fz
        #        self.ea = libs.ea
        #
        #        self.temperature = libs.temperature
        #        self.pressure = libs.pressure
        #        self.ekinetic = libs.ekinetic
        #        self.epotential = libs.energy

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

        libs.atom = self.sys.atom
        libs.natom = self.sys.natom
        libs.cell = self.sys.cell
        libs.volume = self.sys.volume
        libs.sites = self.sys.sites
        libs.nsites = self.sys.nsites
        libs.nmolecules = len(self.sys.molecule)
        libs.molecules = self.sys.molecule

#        libs.mass = self.sys.mat
#        libs.rx = self.sys.rx
#        libs.ry = self.sys.ry
#        libs.rz = self.sys.rz
#        libs.vx = self.sys.vx
#        libs.vy = self.sys.vy
#        libs.vz = self.sys.vz
#        libs.fx = self.sys.fx
#        libs.fy = self.sys.fy
#        libs.fz = self.sys.fz
#        libs.ea = self.sys.ea
#        libs.atp = self.sys.atp

        libs.prepare()
        libs.npt_berendsen()

    def hook(self, system):

        self.sys = system
        self.nfree = 3*(system.natom-1)

# retirar esta parte!
##################################
        self.mat = self.sys.mat
        self.atp = self.sys.atp 
        self.chg = self.sys.chg
        self.rx = self.sys.rx
        self.ry = self.sys.ry
        self.rz = self.sys.rz
        self.vx = self.sys.vx
        self.vy = self.sys.vy
        self.vz = self.sys.vz
        self.fx = self.sys.fx
        self.fy = self.sys.fy
        self.fz = self.sys.fz
        self.ea = self.sys.ea
#################################3

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
        
# retirar esta parte!
##################################
        self.sys.mat = self.mat
        self.sys.atp = self.atp 
        self.sys.chg = self.chg
        self.sys.rx = self.rx
        self.sys.ry = self.ry
        self.sys.rz = self.rz
        self.sys.vx = self.vx
        self.sys.vy = self.vy
        self.sys.vz = self.vz
        self.sys.fx = self.fx
        self.sys.fy = self.fy
        self.sys.fz = self.fz
        self.sys.ea = self.ea
##################################

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

