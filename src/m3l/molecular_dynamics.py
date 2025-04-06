import numpy as np
from m3l.structure import Atom, System2
from m3l.utils import Constants
from m3l import libs

class ForceField(Constants):

    params = np.array([])

    def van_der_waals(self, vars_):

        prms = []
        for list_ in vars_:
            prms.append(list_)

        self.params = np.array([prms])

    def interaction(self, system, params):

        libs.rx = system.rx
        libs.ry = system.ry
        libs.rz = system.rz
        libs.natom = len(system.zat)
        libs.cell = system.cell
        libs.params = params
        libs.forces()

        system.epotential = libs.energy

        self.fx = libs.fx 
        self.fy = libs.fy 
        self.fz = libs.fz 
        self.ea = libs.ea 

        return system

    def __call__(self):

        return self.params

class Ensemble(Constants):

    def __init__(self, temp_bath, press_bath, timestep, force_field, tstat = None, pstat = None, bfactor = None):

        self.timestep = np.array(timestep*self.TIMECONV, dtype = np.float32)
        self.dtimestep = np.array(timestep*self.TIMECONV, dtype = np.float32)
        self.force_field = np.array(force_field, dtype = np.float32)
        self.temp_bath = np.array(temp_bath*self.TEMPCONV, dtype = np.float32)
        self.press_bath = np.array(press_bath*self.PCONV, dtype = np.float32)
        self.tstat = np.array(tstat*self.TIMECONV, dtype = np.float32)
        self.pstat = np.array(pstat*self.TIMECONV, dtype = np.float32)
        self.friction = 0.e0

        if bfactor:
            self.bfactor = np.array(bfactor, dtype = np.float32)
        else:
            self.bfactor = self.BETAFACTOR
    
        self.bfactor = self.bfactor*self.PCONVINV

    def nve_verlet(self):

        print("não valendo...")
        #        ForceField().interaction(self.system, self.timestep, self.force_field)
        #
        #        for atom in self.system.atoms:
        #            mass = Atom().setMass(atom[0])
        #            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
        #            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
        #            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass
        #            atom[1] = atom[1]+self.timestep*atom[4]
        #            atom[2] = atom[2]+self.timestep*atom[5]
        #            atom[3] = atom[3]+self.timestep*atom[6]
        #
        #        self.ccp()
        #
        #        ForceField().interaction(self.system, self.timestep, self.force_field)
        #
        #        for atom in self.system.atoms:
        #            mass = Atom().setMass(atom[0])
        #            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
        #            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
        #            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass

    def nvt_verlet(self):

        print("não valendo...")

    def npt_verlet(self):

        libs.sigma = 0.5*self.nfree*self.KB*self.temp_bath
        libs.natom = len(self.mat)
        libs.timestep = self.dtimestep
        libs.bfactor = self.bfactor
        libs.tstat = self.tstat
        libs.pstat = self.pstat
        libs.press_bath = self.press_bath
        libs.cell = self.cell

        libs.allocate_arrays()

        libs.params = self.force_field 

        libs.mass = self.mat
        libs.rx = self.rx
        libs.ry = self.ry
        libs.rz = self.rz
        libs.vx = self.vx
        libs.vy = self.vy
        libs.vz = self.vz
        libs.atp = self.atp

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

#        self.setTimestep(libs.drmax)
#
#    def setTimestep(self, drmax):
#
#        if drmax >= 0.5*self.ACONV:
#            self.dtimestep = 0.95*self.dtimestep
#        else:
#            self.dtimestep = self.timestep

    def hook(self, system):

#        self.cell = system.cell
#        self.temperature = system.temperature
#        self.pressure = system.pressure
#        self.epotential = system.epotential
#        self.ekinetic = system.ekinetic
#        self.nfree = 3*(len(system.mat)-1)
#
#        self.mat = system.mat
#        self.atp = system.atp
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

        self.cell = np.array(system.cell)
        self.temperature = np.array(system.temperature)
        self.pressure = np.array(system.pressure)
        self.epotential = np.array(system.epotential)
        self.ekinetic = np.array(system.ekinetic)
        self.nfree = 3*(len(system.mat)-1)

        self.mat = np.array(system.mat)
        self.atp = np.array(system.atp)
        self.rx = np.array(system.rx)
        self.ry = np.array(system.ry)
        self.rz = np.array(system.rz)
        self.vx = np.array(system.vx)
        self.vy = np.array(system.vy)
        self.vz = np.array(system.vz)
        self.fx = np.array(system.fx)
        self.fy = np.array(system.fy)
        self.fz = np.array(system.fz)
        self.ea = np.array(system.ea)

    def hook_output(self):

        context = System2(
                self.cell,
                self.temperature,
                self.pressure,
                self.epotential,
                self.ekinetic,
                self.mat,
                self.atp,
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
        if self.tstat:
            if self.pstat:
                self.npt_verlet()
            else:
                self.nvt_verlet()
        else:
            self.nve_verlet()
        output = self.hook_output()
        return output 

