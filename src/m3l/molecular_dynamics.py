from numpy import sqrt
from m3l.structure import System 

class ForceField():

    def __init__(self):

        self.k_cte = 1.25
        self.r_cte = 1.88

    def twoPaticle(self, atoms):

        epot = 0.0
        for i in range(len(atoms)):
            for j in range(i+1,len(atoms)):
                x = atoms[j]['x']-atoms[i]['x']
                y = atoms[j]['y']-atoms[i]['y']
                z = atoms[j]['z']-atoms[i]['z']
                dr = sqrt(x**2+y**2+z**2)
                epot += 0.5*self.k_cte*(dr-self.r_cte)**2
                fr = self.k_cte*(dr-self.r_cte)/dr
                atoms[i]['fx']+=fr*x
                atoms[i]['fy']+=fr*y
                atoms[i]['fz']+=fr*z
                atoms[j]['fx']+=-fr*x
                atoms[j]['fy']+=-fr*y
                atoms[j]['fz']+=-fr*z

        return atoms

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
            self.atoms = ForceField().twoPaticle(self.atoms)
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
                'energy': self.kinetic_energy*self.ECONV,
                'temperature': self.temperature*self.TEMPCONV,
                's2': 0.e0
                })

    def setVelocity(self):

        nfree = 3*(self.natoms-1)

        for atom in self.atoms:
            atom['vx'] = sqrt(nfree*self.KB*self.temperature/(3.0*atom['mass']))
            atom['vy'] = sqrt(nfree*self.KB*self.temperature/(3.0*atom['mass']))
            atom['vz'] = sqrt(nfree*self.KB*self.temperature/(3.0*atom['mass']))

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
