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
        libs.natom = self.natom
        libs.timestep = self.dtimestep
        libs.bfactor = self.bfactor
        libs.tstat = self.tstat
        libs.pstat = self.pstat
        libs.press_bath = self.press_bath
        libs.temp_bath = self.temp_bath
        libs.params = self.force_field 
        libs.nvdw = len(self.force_field)-1
        libs.rvdw = self.force_field[0, 0]

        libs.cell = self.sys.cell
        libs.volume = self.sys.volume
        libs.sites = self.sys.sites
        libs.nsites = self.sys.nsites
        libs.nmolecules = len(self.sys.molecules)
        libs.molecules = self.sys.molecules

        libs.mass = self.sys.mat
        libs.rx = self.sys.rx
        libs.ry = self.sys.ry
        libs.rz = self.sys.rz
        libs.vx = self.sys.vx
        libs.vy = self.sys.vy
        libs.vz = self.sys.vz
        libs.fx = self.sys.fx
        libs.fy = self.sys.fy
        libs.fz = self.sys.fz
        libs.ea = self.sys.ea
        libs.atp = self.sys.atp

        libs.prepare()
        libs.npt_berendsen()

        self.updateSystem()

    def updateSystem(self):

        description = self.sys.description

        self.sys = System()

        self.sys.description = description
        self.sys.cell = libs.cell
        self.sys.volume = libs.volume
        self.sys.sites = libs.sites
        self.sys.nsites = libs.nsites
        self.sys.molecules = libs.molecules
        
        self.sys.mat = libs.mass
        self.sys.atp = libs.atp 
        self.sys.chg = self.chg
        self.sys.rx = libs.rx
        self.sys.ry = libs.ry
        self.sys.rz = libs.rz
        self.sys.vx = libs.vx
        self.sys.vy = libs.vy
        self.sys.vz = libs.vz
        self.sys.fx = libs.fx
        self.sys.fy = libs.fy
        self.sys.fz = libs.fz
        self.sys.ea = libs.ea

        self.sys.temperature = libs.temperature
        self.sys.temp_friction = libs.temp_friction
        self.sys.temp_bath = libs.temp_bath
        self.sys.pressure = libs.pressure
        self.sys.press_bath = libs.press_bath
        self.sys.press_friction = libs.press_friction
        self.sys.ekinetic = libs.ekinetic
        self.sys.epotential = libs.energy
        self.sys.volume = libs.volume

    def hook(self, system):

        self.sys = system
#        self.description = system.description
#        self.cell = system.cell
#        self.volume = system.volume
#        self.temperature = system.temperature
#        self.temp_friction = system.temp_friction
#        self.pressure = system.pressure
#        self.press_friction = system.press_friction
#        self.epotential = system.epotential
#        self.ekinetic = system.ekinetic
        self.natom = len(self.sys.mat)
        self.nfree = 3*(self.natom-1)
#        self.nsites = system.nsites
#        self.sites = system.sites
#        self.molecules = system.molecules
#
#        self.mat = system.mat
#        self.atp = system.atp
        self.chg = system.chg
#        self.rx = system.rx
#        self.ry = system.ry
#        self.rz = system.rz
#        self.vx = system.vx
#        self.vy = system.vy
#        self.vz = system.vz
#        self.fx = system.fx
#        self.fy = system.fy
#        self.fz = system.fz
#        self.ea = system.ea

    def hook_output(self):

        #        context = System2(
        #                self.description,
        #                self.cell,
        #                self.volume,
        #                self.temperature,
        #                self.temp_friction,
        #                self.pressure,
        #                self.press_friction,
        #                self.epotential,
        #                self.ekinetic,
        #                self.nsites,
        #                self.sites,
        #                self.molecules,
        #                self.mat,
        #                self.atp,
        #                self.chg,
        #                self.rx,
        #                self.ry,
        #                self.rz,
        #                self.vx,
        #                self.vy,
        #                self.vz,
        #                self.fx,
        #                self.fy,
        #                self.fz,
        #                self.ea,
        #                )

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

