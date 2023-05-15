from numpy import sqrt
from m3l.structure import System 

class ForceField():

    def __init__(self):

        super().__init__()

        return

class Integration(System):

#    def __init__(self):

#        super().__init__(acell, bcell, ccell, filename)

    timestep = 0.0

    def nve(self):

        for atom in self.atoms:
            atom['vx'] = atom['vx']+0.5*self.timestep/atom['mass']
            atom['vy'] = atom['vy']+0.5*self.timestep/atom['mass']
            atom['vz'] = atom['vz']+0.5*self.timestep/atom['mass']
            atom['x'] = atom['x']+self.timestep*atom['vx']
            atom['y'] = atom['y']+self.timestep*atom['vy']
            atom['z'] = atom['z']+self.timestep*atom['vz']

        self.ccp()

        for atom in self.atoms:
            atom['vx'] = atom['vx']+0.5*self.timestep/atom['mass']
            atom['vy'] = atom['vy']+0.5*self.timestep/atom['mass']
            atom['vz'] = atom['vz']+0.5*self.timestep/atom['mass']

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

    def convertUnitsInv(self):

        self.acell = self.acell*self.aconv
        self.bcell = self.bcell*self.aconv
        self.ccell = self.ccell*self.aconv
        self.timestep = self.timestep*self.timeconv

        for at in self.atoms:
            at['x'] = at['x']*self.aconv
            at['y'] = at['y']*self.aconv
            at['z'] = at['z']*self.aconv
            at['mass'] = at['mass']*self.mconv

    def convertUnits(self):

        self.acell = self.acell/self.aconv
        self.bcell = self.bcell/self.aconv
        self.ccell = self.ccell/self.aconv
        self.timestep = self.timestep/self.timeconv

        for at in self.atoms:
            at['x'] = at['x']/self.aconv
            at['y'] = at['y']/self.aconv
            at['z'] = at['z']/self.aconv
            at['mass'] = at['mass']/self.mconv

    def running(self, timestep, nstep, temperature):

        self.timestep = timestep
        self.nsteps = nstep
        self.temperature_input = temperature

        self.setVelocity()

        for step in range(1,nstep+1):
            self.convertUnits()
            self.ccp()
            self.nve()
            self.setTemperature()
            self.convertUnitsInv()
            self.setFrame(step)

        return

    def setFrame(self, step):

        for at in self.atoms:
            self.frames.append({
                'step': step,
                'id': at['id'],
                'mass': 1.0,
                'x': at['x'],
                'y': at['y'],
                'z': at['z'],
                'energy': 0.e0,
                'temperature': self.temperature,
                's2': 0.e0
                })

    def setVelocity(self):

        for at in self.atoms:
            at['vx'] = sqrt(self.kb*self.temperature_input/at['mass'])
            at['vy'] = sqrt(self.kb*self.temperature_input/at['mass'])
            at['vz'] = sqrt(self.kb*self.temperature_input/at['mass'])

    def setTemperature(self):

        nfree = 3*(self.natoms-1)

        soma = 0.0
        for atom in self.atoms:
            soma += atom['mass']*(atom['vx']**2+atom['vy']**2+atom['vz']**2)

        self.temperature = soma/(self.kb*nfree)

    def __str__(self):
        return (f"""---Molecular dynamics module---
                System:
                - Orthorrombic cell: {self.volume['value']*self.aconv**3} {self.volume['unit']}
                - Total of atoms: {self.natoms}""")
