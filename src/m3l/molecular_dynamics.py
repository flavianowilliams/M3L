import numpy as np
from m3l.structure import Atom, System2
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
            self.tstat = np.array(0.1*self.TIMECONV, dtype = np.float64)

        if pstat:
            self.pstat = np.array(pstat*self.TIMECONV, dtype = np.float64)
        else:
            self.pstat = np.array(0.1*self.TIMECONV, dtype = np.float64)

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

#        libs.sigma = 0.5*self.nfree*self.KB*self.temp_bath
        libs.nfree = self.nfree
        libs.temp_bath = self.temp_bath
        libs.natom = len(self.mat)
        libs.timestep = self.dtimestep
        libs.tstat = self.tstat
        libs.cell = self.cell
        libs.sites = self.sites
        libs.nsites = self.nsites
        libs.nmolecules = len(self.molecules)
        libs.molecules = self.molecules

        libs.params = self.force_field 
        libs.nvdw = len(self.force_field)-1
        libs.rvdw = self.force_field[0, 0]

        libs.mass = self.mat
        libs.rx = self.rx
        libs.ry = self.ry
        libs.rz = self.rz
        libs.vx = self.vx
        libs.vy = self.vy
        libs.vz = self.vz
        libs.fx = self.fx
        libs.fy = self.fy
        libs.fz = self.fz
        libs.ea = self.ea
        libs.atp = self.atp

        libs.prepare()
        libs.nvt_berendsen()

        self.cell = libs.cell
        
        self.rx = libs.rx
        self.ry = libs.ry
        self.rz = libs.rz
        self.vx = libs.vx
        self.vy = libs.vy
        self.vz = libs.vz
        self.fx = libs.fx
        self.fy = libs.fy
        self.fz = libs.fz
        self.ea = libs.ea

        self.temperature = libs.temperature
        self.pressure = libs.pressure
        self.ekinetic = libs.ekinetic
        self.epotential = libs.energy

    def npt_verlet(self):

#        libs.sigma = 0.5*self.nfree*self.KB*self.temp_bath
        libs.nfree = self.nfree
        libs.natom = len(self.mat)
        libs.timestep = self.dtimestep
        libs.bfactor = self.bfactor
        libs.tstat = self.tstat
        libs.pstat = self.pstat
        libs.cell = self.cell
        libs.sites = self.sites
        libs.nsites = self.nsites
        libs.nmolecules = len(self.molecules)
        libs.molecules = self.molecules
        libs.press_bath = self.press_bath
        libs.temp_bath = self.temp_bath

        libs.params = self.force_field 
        libs.nvdw = len(self.force_field)-1
        libs.rvdw = self.force_field[0, 0]

        libs.mass = self.mat
        libs.rx = self.rx
        libs.ry = self.ry
        libs.rz = self.rz
        libs.vx = self.vx
        libs.vy = self.vy
        libs.vz = self.vz
        libs.fx = self.fx
        libs.fy = self.fy
        libs.fz = self.fz
        libs.ea = self.ea
        libs.atp = self.atp

        libs.prepare()
        libs.npt_berendsen()

        self.cell = libs.cell
        
        self.rx = libs.rx
        self.ry = libs.ry
        self.rz = libs.rz
        self.vx = libs.vx
        self.vy = libs.vy
        self.vz = libs.vz
        self.fx = libs.fx
        self.fy = libs.fy
        self.fz = libs.fz
        self.ea = libs.ea

        self.temperature = libs.temperature
        self.pressure = libs.pressure
        self.ekinetic = libs.ekinetic
        self.epotential = libs.energy

    def hook(self, system):

        self.description = system.description
        self.cell = system.cell
        self.temperature = system.temperature
        self.pressure = system.pressure
        self.epotential = system.epotential
        self.ekinetic = system.ekinetic
        self.nfree = 3*(len(system.mat)-1)
        self.nsites = system.nsites
        self.sites = system.sites
        self.molecules = system.molecules

        self.mat = system.mat
        self.atp = system.atp
        self.chg = system.chg
        self.rx = system.rx
        self.ry = system.ry
        self.rz = system.rz
        self.vx = system.vx
        self.vy = system.vy
        self.vz = system.vz
        self.fx = system.fx
        self.fy = system.fy
        self.fz = system.fz
        self.ea = system.ea

    def hook_output(self):

        context = System2(
                self.description,
                self.cell,
                self.temperature,
                self.pressure,
                self.epotential,
                self.ekinetic,
                self.nsites,
                self.sites,
                self.molecules,
                self.mat,
                self.atp,
                self.chg,
                self.rx,
                self.ry,
                self.rz,
                self.vx,
                self.vy,
                self.vz,
                self.fx,
                self.fy,
                self.fz,
                self.ea,
                )

        return context

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

