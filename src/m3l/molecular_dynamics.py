from numpy import sqrt
from m3l.structure import System
from m3l.utils import Conversion 

class ForceField(Conversion):

    def __init__(self, atoms):

        self.k_cte = 15/(self.ECONV/self.ACONV**2)
        self.r_cte = 1.0/self.ACONV
        self.atm = atoms

        self.twoPaticle()

    def twoPaticle(self):

        self.epot = 0.0
        for i in range(len(self.atm)):
            for j in range(i+1,len(self.atm)):
                x = self.atm[j]['x']-self.atm[i]['x']
                y = self.atm[j]['y']-self.atm[i]['y']
                z = self.atm[j]['z']-self.atm[i]['z']
                dr = sqrt(x**2+y**2+z**2)
                self.epot += 0.5*self.k_cte*(dr-self.r_cte)**2
                fr = -self.k_cte*(dr-self.r_cte)/dr
                self.atm[i]['fx']+=-fr*x
                self.atm[i]['fy']+=-fr*y
                self.atm[i]['fz']+=-fr*z
                self.atm[j]['fx']+=fr*x
                self.atm[j]['fy']+=fr*y
                self.atm[j]['fz']+=fr*z

    def getAtm(self):
        return self.atm

    def getEpot(self):
        return self.epot 

class Integration(System):

#    def __init__(self):

#        super().__init__(acell, bcell, ccell, filename)

    timestep = 0.0

    def nve(self):

        for atom in self.atoms:
            atom['vx'] = atom['vx']+0.5*self.timestep*atom['fx']/atom['mass']
            atom['vy'] = atom['vy']+0.5*self.timestep*atom['fy']/atom['mass']
            atom['vz'] = atom['vz']+0.5*self.timestep*atom['fz']/atom['mass']
            atom['x'] = atom['x']+self.timestep*atom['vx']
            atom['y'] = atom['y']+self.timestep*atom['vy']
            atom['z'] = atom['z']+self.timestep*atom['vz']

        self.ccp()

        for atom in self.atoms:
            atom['vx'] = atom['vx']+0.5*self.timestep*atom['fx']/atom['mass']
            atom['vy'] = atom['vy']+0.5*self.timestep*atom['fy']/atom['mass']
            atom['vz'] = atom['vz']+0.5*self.timestep*atom['fz']/atom['mass']

        return

class MolecularDynamics(Integration):

    frames = []

    def __init__(self, acell, bcell, ccell, filename):

        super().__init__(acell, bcell, ccell, filename)

        self.setAtoms()
        self.setNAtoms()
        self.setVolume()
        self.ccp()

        return

    def convertUnits(self):

        self.acell = self.acell/self.ACONV
        self.bcell = self.bcell/self.ACONV
        self.ccell = self.ccell/self.ACONV
        self.timestep = self.timestep/self.TIMECONV
        self.temperature = self.temperature/self.TEMPCONV

        for atom in self.atoms:
            atom['x'] = atom['x']/self.ACONV
            atom['y'] = atom['y']/self.ACONV
            atom['z'] = atom['z']/self.ACONV
            atom['mass'] = atom['mass']/self.MCONV

    def running(self, timestep, nstep, temperature):

        self.timestep = timestep
        self.nsteps = nstep
        self.temperature = temperature

        self.convertUnits()

        self.setVelocity()

        for step in range(1,nstep+1):
            self.ccp()
            force_field = ForceField(self.atoms)
            self.atoms = force_field.getAtm()
            self.epot = force_field.getEpot()
            self.nve()
            self.setKineticEnergy()
            self.setTemperature()
            self.setFrame(step)

        return

    def setFrame(self, step):

        for at in self.atoms:
            self.frames.append({
                'step': step,
                'id': at['id'],
                'mass': at['mass']*self.MCONV,
                'x': at['x']*self.ACONV,
                'y': at['y']*self.ACONV,
                'z': at['z']*self.ACONV,
                'energy': (self.kinetic_energy+self.epot)*self.ECONV,
                'temperature': self.temperature*self.TEMPCONV,
                's2': 0.e0
                })

    def setVelocity(self):

        nfree = 3*(self.natoms-1)

        for atom in self.atoms:
            atom['vx'] = sqrt(nfree*self.KB*self.temperature/(6.0*atom['mass']))
            atom['vy'] = sqrt(nfree*self.KB*self.temperature/(6.0*atom['mass']))
            atom['vz'] = sqrt(nfree*self.KB*self.temperature/(6.0*atom['mass']))

    def setKineticEnergy(self):

        soma = 0.0
        for atom in self.atoms:
            soma += atom['mass']*(atom['vx']**2+atom['vy']**2+atom['vz']**2)

        self.kinetic_energy = 0.5*soma

    def setTemperature(self):

        nfree = 3*(self.natoms-1)

        self.temperature = 2.0*self.kinetic_energy/(self.KB*nfree)

    def __str__(self):
        return (f"""---Molecular dynamics module---
                System:
                - Orthorrombic cell: {self.volume['value']} {self.volume['unit']}
                - Total of atoms: {self.natoms}""")
