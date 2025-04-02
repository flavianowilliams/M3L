from m3l.molecular_dynamics import ForceField
from m3l.structure import Constants
from m3l import libs
import numpy as np

class Optimize(Constants):

    def __init__(self, force_field, learning_rate):

        self.force_field = force_field
        self.learning_rate = learning_rate

    def SD(self, system):

        ForceField().interaction(system, self.force_field)

        for atom in system.atoms:
            atom[1] = atom[1] + self.learning_rate*atom[7]
            atom[2] = atom[2] + self.learning_rate*atom[8]
            atom[3] = atom[3] + self.learning_rate*atom[9]

        rx = []
        ry = []
        rz = []
        for atom in system.atoms:
            rx.append(atom[1])
            ry.append(atom[2])
            rz.append(atom[3])

        libs.natom = len(system.atoms)
        libs.cell = system.cell
        libs.rx = rx
        libs.ry = ry
        libs.rz = rz

        libs.ccp()

        for i, atom in enumerate(system.atoms):
            atom[1] = np.array(libs.rx[i])
            atom[2] = np.array(libs.ry[i])
            atom[3] = np.array(libs.rz[i])

        return system

