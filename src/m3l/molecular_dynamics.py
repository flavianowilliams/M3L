import numpy as np
from m3l.structure import Atom
from m3l.utils import Constants
from m3l import libs

class ForceField(Constants):

    params = []
    r_cutoff = np.array(0.0)

    def interaction(self, system, timestep, params):

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
        libs.params = params
        libs.cell = system.cell
        libs.timestep = timestep
        libs.forces()

        for i, atom in enumerate(system.atoms):
            atom[7] = np.array(libs.fx[i])
            atom[8] = np.array(libs.fy[i])
            atom[9] = np.array(libs.fz[i])
            atom[10] = np.array(libs.ea[i])

        return system.atoms

    def __call__(self):

        return np.array(self.params)

class Thermodynamics(Constants):

    def __init__(self, temp_ext, press_ext):
        
        self.temp_ext = temp_ext
        self.press_ext = press_ext

    def setKineticEnergy(self, atoms):

        soma = 0.0
        for atom in atoms:
            mass = Atom().setMass(atom[0])
            soma += mass*(atom[4]**2+atom[5]**2+atom[6]**2)

        self.kinetic_energy = np.array(0.5*soma)

        return self.kinetic_energy

    def setTemperature(self, atoms):

        nfree = 3*(len(atoms)-1)

        self.temperature = 2.0*self.kinetic_energy/(self.KB*nfree)
        self.temperature = np.array(self.temperature)

        return self.temperature 

class Ensemble(Thermodynamics):

    def __init__(self, temp_ext, press_ext, timestep, force_field, integrator):
        super().__init__(temp_ext, press_ext)

        self.timestep = timestep
        self.temp_scale = np.array(0.2)/self.TIMECONV
        self.integrator = integrator
        self.force_field = force_field
    
    def nve_verlet(self):

        ForceField().interaction(self.system, self.timestep, self.force_field)

        for atom in self.system.atoms:
            mass = Atom().setMass(atom[0])
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass
            atom[1] = atom[1]+self.timestep*atom[4]
            atom[2] = atom[2]+self.timestep*atom[5]
            atom[3] = atom[3]+self.timestep*atom[6]

        self.ccp()

        ForceField().interaction(self.system, self.timestep, self.force_field)

        for atom in self.system.atoms:
            mass = Atom().setMass(atom[0])
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass

        return self.system

    def nvt_verlet(self):

        self.sigma = 0.5*self.nfree*self.KB*self.temp_ext

        ForceField().interaction(self.system, self.timestep, self.force_field)

        for atom in self.system.atoms:
            mass = Atom().setMass(atom[0])
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass
            atom[1] = atom[1]+self.timestep*atom[4]
            atom[2] = atom[2]+self.timestep*atom[5]
            atom[3] = atom[3]+self.timestep*atom[6]

        self.ccp()

        ForceField().interaction(self.system, self.timestep, self.force_field)

        for atom in self.system.atoms:
            mass = Atom().setMass(atom[0])
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass

        self.setKineticEnergy(self.system.atoms)

        for atom in self.system.atoms:
            atom[4] = atom[4]*self.scaling()
            atom[5] = atom[5]*self.scaling()
            atom[6] = atom[6]*self.scaling()

        self.setKineticEnergy(self.system.atoms)
        self.setTemperature(self.system.atoms)

        self.system.temperature = self.temperature

        return self.system

    def scaling(self):
        value = np.sqrt(1.0+self.timestep*(self.sigma/self.kinetic_energy-1.0)/self.temp_scale)
        return value

    def ccp(self):
        for at in self.system.atoms:
            at[1]=at[1]-self.system.cell[0]*int(at[1]/self.system.cell[0])
            at[2]=at[2]-self.system.cell[1]*int(at[2]/self.system.cell[1])
            at[3]=at[3]-self.system.cell[2]*int(at[3]/self.system.cell[2])
#        else:
#            logging.critical('Constant lattice must be larger than zero!')
#            sys.exit()

    def __call__(self, system):
        self.system = system
        self.nfree = 3*(len(system.atoms)-1)
        if self.integrator == 'nve_verlet':
            return self.nve_verlet()
        elif self.integrator == 'nvt_verlet':
            return self.nvt_verlet()

