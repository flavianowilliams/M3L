from m3l.structure import System as sys
import json
#
system = sys()
system.loadSystem('system_Ar.json')
system.convertUnits()

from m3l.molecular_dynamics import ForceField
class Forces(ForceField):
    def __init__(self):
        super().__init__()
        self.parameters((0.2385/self.ECONV, 3.405/self.ACONV, 10.0/self.ACONV))

model = Forces()

from m3l.optimize import Optimize as optim
#
opt = optim(model(), learning_rate = 1.e-3)

import time
import csv
#
n_steps = 1
start = time.time()
with open('optimize.csv', 'w', newline = '') as csvfile:
    fieldnames = ['step', 'energy']
    writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
    writer.writeheader()
    for step in range(1, n_steps+1):
        system = opt.SD(system)
        writer.writerow({
            'step': step,
            'energy': system.epotential*system.ECONV})
#
end = time.time()
print(f'Elapsed time: {end - start}')

