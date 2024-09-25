from m3l.molecular_dynamics import ForceField
from m3l.structure import Constants

class Optimize(Constants):

    def __init__(self, force_field, learning_rate):

        self.force_field = force_field
        self.learning_rate = learning_rate

    def SD(self, system):

        system.temperature = 0.0
        system.friction = 0.0

        for atom in system.atoms:
            atom[4] = 0.0
            atom[5] = 0.0
            atom[6] = 0.0

        ForceField().interaction(system, self.force_field)

        for atom in system.atoms:
            atom[1] = atom[1] + self.learning_rate*atom[7]
            atom[2] = atom[2] + self.learning_rate*atom[8]
            atom[3] = atom[3] + self.learning_rate*atom[9]

        self.ccp(system)

        return system

    def ccp(self, system):
        for atom in system.atoms:
            atom[1]=atom[1]-system.cell[0]*int(atom[1]/system.cell[0]) 
            atom[2]=atom[2]-system.cell[1]*int(atom[2]/system.cell[1])
            atom[3]=atom[3]-system.cell[2]*int(atom[3]/system.cell[2]) 

