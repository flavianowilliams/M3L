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
        self.rx_new = []
        self.ry_new = []
        self.rz_new = []
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

            self.rx = self.rx_new
            self.ry = self.ry_new
            self.rz = self.rz_new
#            for i, atom in enumerate(self.atoms):
#                atom[1] = self.rx[i]
#                atom[2] = self.ry[i]
#                atom[3] = self.rz[i]

            self.epotential += de
            self.naccpt += 1
            
        elif np.exp(-(deltvb+self.press_ext*deltvol)) > rnd:

            for i in range(3):
                self.cell[i]=self.celln[i]

            self.rx = self.rx_new
            self.ry = self.ry_new
            self.rz = self.rz_new
#            for i, atom in enumerate(self.atoms):
#                atom[1] = self.rx[i]
#                atom[2] = self.ry[i]
#                atom[3] = self.rz[i]

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

        self.rx_new.clear()
        self.ry_new.clear()
        self.rz_new.clear()

        for i in range(len(self.rx)):
            rnd = np.random.random_sample()
            self.rx_new.append(self.rx[i]+self.drmax*(2.0*rnd-1.0))
            rnd = np.random.random_sample()
            self.ry_new.append(self.ry[i]+self.drmax*(2.0*rnd-1.0))
            rnd = np.random.random_sample()
            self.rz_new.append(self.rz[i]+self.drmax*(2.0*rnd-1.0))

    def energy(self):

        self.setCoordinates()

        libs.natom = len(self.mat)
        libs.cell = self.cell
        libs.rx = self.rx_new
        libs.ry = self.ry_new
        libs.rz = self.rz_new
        libs.params = self.force_field 

        libs.forces()

        return (libs.energy - self.epotential)

    def hook(self, system):

        #        atm = []
        #        for atom in system.atoms:
        #            atm.append(atom)
       
#        self.atoms = np.array(atm)
        self.cell = np.array(system.cell)
        self.temperature = np.array(system.temperature)
        self.pressure = np.array(system.pressure)
        self.epotential = np.array(system.epotential)
        self.ekinetic = np.array(system.ekinetic)
        self.nfree = 3*(len(system.mat)-1)
        self.mat = np.array(system.mat)
        self.atp = np.array(system.atp)
        self.rx = np.array(system.rx)
        self.ry = np.array(system.ry)
        self.rz = np.array(system.rz)
        self.vx = np.array(system.vx)
        self.vy = np.array(system.vy)
        self.vz = np.array(system.vz)
        self.fx = np.array(system.fx)
        self.fy = np.array(system.fy)
        self.fz = np.array(system.fz)
        self.ea = np.array(system.ea)


    def hook_output(self):

        context = System2(
                self.cell,
                self.temperature,
                self.pressure,
                self.epotential,
                self.ekinetic,
                self.mat,
                self.atp,
                self.rx,
                self.ry,
                self.rz,
                self.vx,
                self.vy,
                self.vz,
                self.fx,
                self.fy,
                self.fz,
                self.ea,
                )

        return context

    def __call__(self, system):
        self.hook(system)
        self.metropolis()
        output = self.hook_output()
        return output
