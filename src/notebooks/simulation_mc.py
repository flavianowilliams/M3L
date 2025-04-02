from m3l.structure import System as sys
import json
#
system = sys()
system.loadSystem('system_new.json')
system.convertUnits()
system.save('system_au.json')
#system.convertUnitsInv()
#system.atoms

# definindo o modelo de interação entre os átomos (campo de força)
from m3l.molecular_dynamics import ForceField
class Forces(ForceField):
    def __init__(self):
        super().__init__()
        self.parameters((125.7/self.ECONV, 3.345/self.ACONV, 9.0/self.ACONV))

model = Forces()

# definindo o modelo estatístico (ensemble)
from m3l.statistics import MonteCarlo as mc
#
temperature = 10.0
temperature = temperature/system.TEMPCONV
pressure = 1.0
drmax=1.0/system.ACONV
dvmax=0.0/system.ACONV**3
force_field = model()
monte_carlo = mc(temperature, drmax, dvmax, force_field, nadjst = 200)
tmp = monte_carlo(system)
#print(system.epotential, tmp.epotential)

# executando looping
import time
import csv
#
n_steps = 100000
i_step = 0
en = []
drmax = []
start = time.time()
with open('history.csv', 'w', newline = '') as csvfile1, open('thermodynamics.csv', 'w', newline = '') as csvfile2:
    fieldnames1 = ['step', 'atom', 'x', 'y', 'z']
    writer1 = csv.DictWriter(csvfile1, fieldnames = fieldnames1)
    writer1.writeheader()
    fieldnames2 = ['step', 'energy']
    writer2 = csv.DictWriter(csvfile2, fieldnames = fieldnames2)
    writer2.writeheader()
    for step in range(n_steps):
        chk1 = monte_carlo.naccpt
        system = monte_carlo(system)
        chk2 = monte_carlo.naccpt
        drmax.append(monte_carlo.drmax)
        if chk2 != chk1:
            volume = system.cell[0]*system.cell[1]*system.cell[2]
            volume = volume*system.ACONV**3
            writer2.writerow({
                'step': step,
                'energy': system.epotential*system.ECONV})
            for atom in system.atoms:
                writer1.writerow({
                    'step': step,
                    'atom': atom[0],
                    'x': atom[1]*system.ACONV,
                    'y': atom[2]*system.ACONV,
                    'z': atom[3]*system.ACONV})
            en.append(system.epotential*system.ECONV)
            i_step += 1
#        print(f'i: {step}; energy: {system.energy*system.ECONV}')
end = time.time()
print(f'Elapsed time: {end - start}')

# salvando sistema no formato JSON
system.convertUnitsInv()
system.save('system_new.json')
