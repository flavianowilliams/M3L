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

#    def convertUnitsInv(self):
#
#       self.acell = self.acell*self.ACONV
#       self.bcell = self.bcell*self.ACONV
#       self.ccell = self.ccell*self.ACONV
#       self.timestep = self.timestep*self.TIMECONV
#       self.temperature = self.temperature*self.TEMPCONV

#       for at in self.atoms:
#           at['x'] = at['x']*self.ACONV
#           at['y'] = at['y']*self.ACONV
#           at['z'] = at['z']*self.ACONV
#           at['mass'] = at['mass']*self.MCONV

    def convertUnits(self):

        self.acell = self.acell/self.ACONV
        self.bcell = self.bcell/self.ACONV
        self.ccell = self.ccell/self.ACONV
        self.timestep = self.timestep/self.TIMECONV
        self.temperature = self.temperature/self.TEMPCONV

        for at in self.atoms:
            at['x'] = at['x']/self.ACONV
            at['y'] = at['y']/self.ACONV
            at['z'] = at['z']/self.ACONV
            at['mass'] = at['mass']/self.MCONV

    def running(self, timestep, nstep, temperature):

        self.timestep = timestep
        self.nsteps = nstep
        self.temperature = temperature

        self.convertUnits()

        self.setVelocity()

        for step in range(1,nstep+1):
            self.ccp()
            self.nve()
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
                'energy': 0.e0,
                'temperature': self.temperature*self.TEMPCONV,
                's2': 0.e0
                })

    def setVelocity(self):

        for at in self.atoms:
            at['vx'] = sqrt(self.KB*self.temperature/at['mass'])
            at['vy'] = sqrt(self.KB*self.temperature/at['mass'])
            at['vz'] = sqrt(self.KB*self.temperature/at['mass'])

    def setTemperature(self):

        nfree = 3*(self.natoms-1)

        soma = 0.0
        for atom in self.atoms:
            soma += atom['mass']*(atom['vx']**2+atom['vy']**2+atom['vz']**2)

        self.temperature = soma/(self.KB*nfree)

    def __str__(self):
        return (f"""---Molecular dynamics module---
                System:
                - Orthorrombic cell: {self.volume['value']} {self.volume['unit']}
                - Total of atoms: {self.natoms}""")
