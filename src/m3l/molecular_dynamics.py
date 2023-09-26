import csv
from numpy import sqrt
from m3l.structure import System
from m3l.utils import Conversion

class ForceField(Conversion):

    def __init__(self, atoms):

        self.e_cte = 120/self.ECONV
        self.sigma_cte = 3.4/self.ACONV
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
                self.epot += 4.0e0*self.e_cte*((dr/self.sigma_cte)**12-(dr/self.sigma_cte)**6)
                fr = 24.0e0*self.e_cte*(2.0e0*(dr/self.sigma_cte)**12-(dr/self.sigma_cte)**6)
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

    def __init__(self, acell, bcell, ccell, filename, **kwargs):

        super().__init__(acell, bcell, ccell, filename)

        if kwargs['reload'] is not None:
            self.reload = kwargs['reload']

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
            atom['vx'] = atom['vx']/(self.ACONV/self.TIMECONV)
            atom['vy'] = atom['vy']/(self.ACONV/self.TIMECONV)
            atom['vz'] = atom['vz']/(self.ACONV/self.TIMECONV)
            atom['fx'] = atom['fx']/(self.ECONV/self.ACONV)
            atom['fy'] = atom['fy']/(self.ECONV/self.ACONV)
            atom['fz'] = atom['fz']/(self.ECONV/self.ACONV)
            atom['mass'] = atom['mass']/self.MCONV

    def running(self, timestep, nstep, temperature):

        self.timestep = timestep
        self.nsteps = nstep
        self.temperature = temperature

        self.convertUnits()

        self.setVelocity()

        with open('history.csv', 'w', newline='') as file:

            fields = ['step', 'id', 'symbol', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'energy', 'temperature', 's2']
            history = csv.DictWriter(file, fieldnames=fields)
            history.writeheader()

            for step in range(1,nstep+1):
                self.ccp()
                force_field = ForceField(self.atoms)
                self.atoms = force_field.getAtm()
                self.epot = force_field.getEpot()
                self.nve()
                self.setKineticEnergy()
                self.setTemperature()
                self.setFrame(step)
                history.writerows(self.frame)

        with open('system.csv', 'w', newline='') as file:

            fields = ['id', 'symbol', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 's2']
            system = csv.DictWriter(file, fieldnames=fields)
            system.writeheader()
            for atom in self.frame:
                system.writerow({
                    'id': atom['id'],
                    'symbol': atom['symbol'],
                    'mass': atom['mass'],
                    'x': atom['x'],
                    'y': atom['y'],
                    'z': atom['z'],
                    'vx': atom['vx'],
                    'vy': atom['vy'],
                    'vz': atom['vz'],
                    'fx': atom['fx'],
                    'fy': atom['fy'],
                    'fz': atom['fz'],
                    's2': atom['s2'],
                    })

        return

    def setFrame(self, step):

        self.frame = []
        for atom in self.atoms:
            self.frame.append({
                'step': step,
                'id': atom['id'],
                'symbol': atom['symbol'],
                'mass': atom['mass']*self.MCONV,
                'x': atom['x']*self.ACONV,
                'y': atom['y']*self.ACONV,
                'z': atom['z']*self.ACONV,
                'vx': atom['vx']*(self.ACONV/self.TIMECONV),
                'vy': atom['vy']*(self.ACONV/self.TIMECONV),
                'vz': atom['vz']*(self.ACONV/self.TIMECONV),
                'fx': atom['fx']*(self.ECONV/self.ACONV),
                'fy': atom['fy']*(self.ECONV/self.ACONV),
                'fz': atom['fz']*(self.ECONV/self.ACONV),
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
