from m3l.structure import System 

class ForceField():

    def __init__(self):

        super().__init__()

        return

class Integration(System):

#    def __init__(self):

#        super().__init__(acell, bcell, ccell, filename)

    def nve(self, step):

        for atom in self.atoms:
            atom['x'] = atom['x']+0.2*self.acell
            atom['y'] = atom['y']+0.0
            atom['z'] = atom['z']+0.0
            atom['step'] = step

        return

class MolecularDynamics(Integration):

    def __init__(self, acell, bcell, ccell, filename):

        super().__init__(acell, bcell, ccell, filename)

        self.setAtoms()
        self.setNAtoms()
        self.convertUnits()
        self.setVolume()
        self.ccp()

        return

    def convertUnitsInv(self):

        self.acell = self.acell*self.aconv
        self.bcell = self.bcell*self.aconv
        self.ccell = self.ccell*self.aconv

        for at in self.atoms:
            at['x'] = at['x']*self.aconv
            at['y'] = at['y']*self.aconv
            at['z'] = at['z']*self.aconv

    def convertUnits(self):

        self.acell = self.acell/self.aconv
        self.bcell = self.bcell/self.aconv
        self.ccell = self.ccell/self.aconv

        for at in self.atoms:
            at['x'] = at['x']/self.aconv
            at['y'] = at['y']/self.aconv
            at['z'] = at['z']/self.aconv

    def running(self, timestep, nstep):

        self.timestep = timestep
        self.nstep = nstep

        for step in range(1,self.nstep+1):
            self.ccp()
            self.nve(step)

        return

    def __str__(self):
        return (f"""---Molecular dynamics module---
                System:
                - Orthorrombic cell: {self.volume['value']*self.aconv**3} {self.volume['unit']}
                - Total of atoms: {self.natoms}""")
