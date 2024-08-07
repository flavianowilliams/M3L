import numpy as np
from m3l.structure import Atom
from m3l.utils import Conversion

class ForceField(Conversion):

    params = []

    def vdw(self, system, params):

        for atom in system.atoms:
            atom[7] = np.array(0.0)
            atom[8] = np.array(0.0)
            atom[9] = np.array(0.0)

        for i in range(len(system.atoms)):
            epot = 0.0
            for j in range(i+1, len(system.atoms)):
                xx = system.atoms[j][1]-system.atoms[i][1]
                yy = system.atoms[j][2]-system.atoms[i][2]
                zz = system.atoms[j][3]-system.atoms[i][3]
                drb = self.bondConstraint(xx, yy, zz, system.cell)
                dr = np.sqrt(drb[0]**2+drb[1]**2+drb[2]**2)
                epot = 4.0*params[0]*((params[1]/dr)**12-(params[1]/dr)**6)
                fr = 24.0*params[0]*(2.0*(params[1]/dr)**12-(params[1]/dr)**6)
                system.atoms[i][7] +=-fr*drb[0]
                system.atoms[i][8] +=-fr*drb[1]
                system.atoms[i][9] +=-fr*drb[2]
                system.atoms[i][10] += epot
                system.atoms[j][7] += fr*drb[0]
                system.atoms[j][8] += fr*drb[1]
                system.atoms[j][9] += fr*drb[2]
                system.atoms[j][10] += epot

        return system.atoms

    def bondConstraint(self, xx, yy, zz, cell):
        xx = xx - cell[0]*int(2.0*(xx/cell[0]))
        yy = yy - cell[1]*int(2.0*(yy/cell[1]))
        zz = zz - cell[2]*int(2.0*(zz/cell[2]))
        return np.array([xx, yy, zz])

    def __call__(self):
        return np.array(self.params)

class Ensemble(Conversion):

    def __init__(self, timestep, force_field, integrator):
        super().__init__()

        self.timestep = timestep
        self.temp_scale = np.array(1.0)/self.TIMECONV
        self.integrator = integrator
        self.force_field = force_field

    def nve(self):

        ForceField().vdw(self.system, self.force_field)

        for atom in self.system.atoms:
            mass = Atom(atom[0]).setMass()
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass
            atom[1] = atom[1]+self.timestep*atom[4]
            atom[2] = atom[2]+self.timestep*atom[5]
            atom[3] = atom[3]+self.timestep*atom[6]

        self.ccp()

        for atom in self.system.atoms:
            mass = Atom(atom[0]).setMass()
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass

        return self.system

    def nvt(self):

        self.force_field.vdw(self.system.atoms)

        for atom in self.system.atoms:
            mass = Atom(atom[0]).setMass()
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass
            atom[1] = atom[1]+self.timestep*atom[4]
            atom[2] = atom[2]+self.timestep*atom[5]
            atom[3] = atom[3]+self.timestep*atom[6]

        self.ccp()

        self.force_field.vdw(self.system.atoms)

        for atom in self.system.atoms:
            mass = Atom(atom[0]).setMass()
            atom[4] = atom[4]+0.5*self.timestep*atom[7]/mass
            atom[5] = atom[5]+0.5*self.timestep*atom[8]/mass
            atom[6] = atom[6]+0.5*self.timestep*atom[9]/mass

    def scaling(self):
        value = np.sqrt(1.0+self.timestep*(self.system.temperature/self.temp_prm-1.0)/self.temp_scale)
        return value

    def ccp(self):
        if self.system.cell[0] > 0.0 and self.system.cell[1] > 0.0 and self.system.cell[2] > 0.0:
            for at in self.system.atoms:
                at[1]=at[1]-self.system.cell[0]*int(at[1]/self.system.cell[0])
                at[2]=at[2]-self.system.cell[1]*int(at[2]/self.system.cell[1])
                at[3]=at[3]-self.system.cell[2]*int(at[3]/self.system.cell[2])
#        else:
#            logging.critical('Constant lattice must be larger than zero!')
#            sys.exit()

    def convertUnits(self):

        self.timestep = self.timestep/self.TIMECONV
        self.temp_prm = self.temp_prm/self.TEMPCONV
#        self.rc_prm = self.rc_prm/self.ACONV

        self.system.cell[0] = self.system.cell[0]/self.ACONV
        self.system.cell[1] = self.system.cell[1]/self.ACONV
        self.system.cell[2] = self.system.cell[2]/self.ACONV

        for atom in self.system.atoms:
            atom[1] = atom[1]/self.ACONV
            atom[2] = atom[2]/self.ACONV
            atom[3] = atom[3]/self.ACONV
            atom[4] = atom[4]/(self.ACONV/self.TIMECONV)
            atom[5] = atom[5]/(self.ACONV/self.TIMECONV)
            atom[6] = atom[6]/(self.ACONV/self.TIMECONV)
            atom[7] = atom[7]/(self.ECONV/self.ACONV)
            atom[8] = atom[8]/(self.ECONV/self.ACONV)
            atom[9] = atom[9]/(self.ECONV/self.ACONV)
            atom[10] = atom[10]/self.ECONV

    def convertUnitsInv(self):

        self.timestep = self.timestep*self.TIMECONV
        self.temp_prm = self.temp_prm*self.TEMPCONV
#        self.rc_prm = self.rc_prm/self.ACONV

        self.system.cell[0] = self.system.cell[0]*self.ACONV
        self.system.cell[1] = self.system.cell[1]*self.ACONV
        self.system.cell[2] = self.system.cell[2]*self.ACONV

        for atom in self.system.atoms:
            atom[1] = atom[1]*self.ACONV
            atom[2] = atom[2]*self.ACONV
            atom[3] = atom[3]*self.ACONV
            atom[4] = atom[4]*(self.ACONV/self.TIMECONV)
            atom[5] = atom[5]*(self.ACONV/self.TIMECONV)
            atom[6] = atom[6]*(self.ACONV/self.TIMECONV)
            atom[7] = atom[7]*(self.ECONV/self.ACONV)
            atom[8] = atom[8]*(self.ECONV/self.ACONV)
            atom[9] = atom[9]*(self.ECONV/self.ACONV)
            atom[10] = atom[10]*self.ECONV

    def __call__(self, system, thermodynamics):
        self.system = system
        self.temp_prm = thermodynamics.temperature
        if self.integrator == 'nve':
            return self.nve()
        elif self.integrator == 'nvt':
            return self.nvt()

class Thermodynamics(Conversion):

    def __init__(self, temperature, pressure):
        super().__init__()

        self.temperature = temperature
        self.pressure = pressure

    def setVelocity(self, atoms):

        for atom in atoms:
            mass = Atom(atom[0]).setMass()
            atom[4] = np.sqrt(self.nfree*self.KB*self.system.temperature/(6.0*mass))
            atom[5] = np.sqrt(self.nfree*self.KB*self.system.temperature/(6.0*mass)) 
            atom[6] = np.sqrt(self.nfree*self.KB*self.system.temperature/(6.0*mass))

    def setKineticEnergy(self, atoms):

        soma = 0.0
        for atom in atoms:
            soma += atom[3]*(atom[4]**2+atom[5]**2+atom[6]**2)

        self.kinetic_energy = np.array(0.5*soma)

        return self.kinetic_energy

    def setTemperature(self):

        temperature = 2.0*self.kinetic_energy/(self.KB*self.nfree)
        temperature = np.array(temperature)

        return self.system.temperature

    def __call__(self, system):
        self.system = system
        self.nfree = (3.0*(len(system.atoms)-1.0))
        self.setKineticEnergy(self.system.atoms)
        self.setTemperature()
        return self.system

#    def history(self):
#
#        with open('history.csv', 'w', newline='') as file:
# 
#            fields = ['id', 'atom', 'time', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'energy', 'temperature', 's2']
#            history = csv.DictWriter(file, fieldnames=fields)
#            history.writeheader()
# 
#            for step in range(1,nstep+1):
#                self.ccp()
#                force_field = ForceField(self.atoms)
##                self.atoms = force_field.getAtm()
#                Symmetry().symm2D()
#                self.epot = force_field.getEpot()
#                Integration(self.timestep).nve()
#                self.setKineticEnergy()
#                self.setTemperature()
#                self.setFrame()
#                history.writerows(self.frame)
#                if step % 50 == 0:
#                    print(f"Step: {step}; Energy: {(self.kinetic_energy+self.epot)*self.ECONV}")
#
#    def __str__(self):
#        return (f"""---Molecular dynamics module---
#                System:
#                - Orthorrombic cell: {self.volume['value']} {self.volume['unit']}
#
#- Total of atoms: {self.natoms}""")
