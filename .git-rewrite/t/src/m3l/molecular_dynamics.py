import numpy as np
from m3l.structure import Atom
from m3l.utils import Constants
from m3l import libs

class ForceField(Constants):

    params = []
    r_cutoff = np.array(0.0)

    def parameters(self, args):

        prms = []
        for list_ in args:
            prms.append(list_)

        self.params = np.array(prms)
    
    def interaction(self, system, params):

        rx = []
        ry = []
        rz = []
        for atom in system.atoms:
            rx.append(atom[1])
            ry.append(atom[2])
            rz.append(atom[3])

        libs.rx = rx
        libs.ry = ry
        libs.rz = rz
        libs.natom = len(system.atoms)
        libs.cell = system.cell
        libs.params = params
        libs.forces()

        system.epotential = libs.energy

        for i, atom in enumerate(system.atoms):
            atom[7] = np.array(libs.fx[i])
            atom[8] = np.array(libs.fy[i])
            atom[9] = np.array(libs.fz[i])
            atom[10] = np.array(libs.ea[i])

        return system

    def __call__(self):

        return self.params

#class Thermodynamics(Constants):
#
#    def __init__(self, temp_ext, press_ext):
#        
#        self.temp_ext = temp_ext
#        self.press_ext = press_ext

class Ensemble(Constants):

    def __init__(self, temp_bath, press_bath, timestep, force_field, tstat = None):

        self.timestep = np.array(timestep)
        self.force_field = force_field
        self.temp_bath = np.array(temp_bath)
        self.press_bath = np.array(press_bath)
        self.tstat = np.array(tstat)
    
#    def setEkinetic(self):
#
#        libs.ekinetic_func()
#
#        return np.array(libs.ekinetic)
#
#    def setTemperature(self):
#
#        ekinetic = self.setEkinetic()
#
#        temperature = 2.0*ekinetic/(self.KB*self.nfree)
#        temperature = np.array(temperature)
#
#        return temperature

    def nve_verlet(self):

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

        return self.system

    def nvt_verlet(self):

        mass = []
        rx = []
        ry = []
        rz = []
        vx = []
        vy = []
        vz = []
        fx = []
        fy = []
        fz = []
        ea = []
        for atom in self.system.atoms:
            mass.append(Atom().setMass(atom[0]))
            rx.append(atom[1])
            ry.append(atom[2])
            rz.append(atom[3])
            vx.append(atom[4])
            vy.append(atom[5])
            vz.append(atom[6])
            fx.append(atom[7])
            fy.append(atom[8])
            fz.append(atom[9])
            ea.append(atom[10])

        libs.sigma = 0.5*self.nfree*self.KB*self.temp_bath
        libs.natom = len(self.system.atoms)
        libs.timestep = self.timestep
        libs.tstat = self.tstat
        libs.friction = self.system.friction
        libs.cell = self.system.cell
        libs.mass = mass
        libs.rx = rx
        libs.ry = ry
        libs.rz = rz
        libs.vx = vx
        libs.vy = vy
        libs.vz = vz
        libs.fx = fx
        libs.fy = fy
        libs.fz = fz
        libs.ea = ea
        libs.params = self.force_field 

        libs.nvt()

        for i, atom in enumerate(self.system.atoms):
            atom[1] = np.array(libs.rx[i])
            atom[2] = np.array(libs.ry[i])
            atom[3] = np.array(libs.rz[i])
            atom[4] = np.array(libs.vx[i])
            atom[5] = np.array(libs.vy[i])
            atom[6] = np.array(libs.vz[i])
            atom[7] = np.array(libs.fx[i])
            atom[8] = np.array(libs.fy[i])
            atom[9] = np.array(libs.fz[i])
            atom[10] = np.array(libs.ea[i])

        self.system.temperature = 2.0*libs.ekinetic/(self.KB*self.nfree)
        self.system.friction = libs.friction
        self.system.ekinetic = libs.ekinetic
        self.system.epotential = libs.energy

        return self.system

#    def ccp(self):
#        for at in self.system.atoms:
#            at[1]=at[1]-self.system.cell[0]*int(at[1]/self.system.cell[0])
#            at[2]=at[2]-self.system.cell[1]*int(at[2]/self.system.cell[1])
#            at[3]=at[3]-self.system.cell[2]*int(at[3]/self.system.cell[2])
#        else:
#            logging.critical('Constant lattice must be larger than zero!')
#            sys.exit()

    def __call__(self, system):
        self.system = system
        self.nfree = 3*(len(system.atoms)-1)
        if self.tstat:
            return self.nvt_verlet()
        else:
            return self.nve_verlet()

