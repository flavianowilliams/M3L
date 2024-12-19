import numpy as np
from m3l import libs
from m3l.structure import System2
from m3l.utils import Constants

class MonteCarlo(Constants):

    ntrial = 0 
    naccpt = 0 
    nadjst = 5

    def __init__(self, temp_ext, press_ext, drmax, dvmax, force_field, **kwargs):

        self.force_field = np.array(force_field)
        self.temp_ext = np.array(temp_ext)
        self.press_ext = np.array(press_ext)
        self.drmax = np.array(drmax)
        self.dvmax = np.array(dvmax)
        self.celln = np.zeros(3)
        self.rx = []
        self.ry = []
        self.rz = []
        self.nadjst = kwargs['nadjst']

    def metropolis(self):

        de = self.energy()
        deltvb = de/(self.KB*self.temp_ext)
        rnd = np.random.random_sample()

        volm = self.cell[0]*self.cell[1]*self.cell[2]
        voln = self.celln[0]*self.celln[1]*self.celln[2]
        deltvol = voln-volm

        self.setDrmax()

        if de < 0.0:

            for i in range(3):
                self.cell[i]=self.celln[i]

            for i, atom in enumerate(self.atoms):
                atom[1] = self.rx[i]
                atom[2] = self.ry[i]
                atom[3] = self.rz[i]

            self.epotential += de
            self.naccpt += 1
            
        elif np.exp(-(deltvb+self.press_ext*deltvol)) > rnd:

            for i in range(3):
                self.cell[i]=self.celln[i]

            for i, atom in enumerate(self.atoms):
                atom[1] = self.rx[i]
                atom[2] = self.ry[i]
                atom[3] = self.rz[i]

            self.epotential += de
            self.naccpt += 1

        self.ntrial += 1

    def setDrmax(self):

        if self.ntrial % self.nadjst == 0:

            ratio = float(self.naccpt)/float(self.nadjst)

            if ratio > 0.5:

                self.drmax = self.drmax*1.05

            else:

                self.drmax = self.drmax*0.95

            self.naccpt = 0

    def setCoordinates(self):

        rnd = np.random.random_sample()
        vn = (self.cell[0]*self.cell[1]*self.cell[2])+(2.0*rnd-1.0)*self.dvmax

        for i in range(3):
            self.celln[i] = vn**(1.0/3.0)

        self.rx.clear()
        self.ry.clear()
        self.rz.clear()

        for atom in self.atoms:
            rnd = np.random.random_sample()
            self.rx.append(atom[1]+self.drmax*(2.0*rnd-1.0))
            rnd = np.random.random_sample()
            self.ry.append(atom[2]+self.drmax*(2.0*rnd-1.0))
            rnd = np.random.random_sample()
            self.rz.append(atom[3]+self.drmax*(2.0*rnd-1.0))

    def energy(self):

        self.setCoordinates()

        libs.natom = len(self.atoms)
        libs.cell = self.cell
        libs.rx = self.rx
        libs.ry = self.ry
        libs.rz = self.rz
        libs.params = self.force_field 

        libs.forces()

        return (libs.energy - self.epotential)

    def hook(self, system):

        atm = []
        for atom in system.atoms:
            atm.append(atom)
       
        self.atoms = np.array(atm)
        self.cell = np.array(system.cell)
        self.temperature = np.array(system.temperature)
        self.friction = np.array(system.friction)
        self.pressure = np.array(system.pressure)
        self.epotential = np.array(system.epotential)
        self.ekinetic = np.array(system.ekinetic)

    def hook_output(self):

        context = System2(
                self.cell,
                self.temperature,
                self.friction,
                self.pressure,
                self.epotential,
                self.ekinetic,
                self.atoms
                )

        return context

    def __call__(self, system):
        self.hook(system)
        self.metropolis()
        output = self.hook_output()
        return output
