import csv
from m3l.symmetry import Symm2D
from m3l.utils import Constant

class DataSet(Constant):

    def __init__(self, eta_prm, rs_prm):

        super().__init__()

        self.eta_prm = eta_prm
        self.rs_prm = rs_prm
        self.atmax = 0
        self.stepmax = 0
        self.atoms = list()
        self.setDS()
        self.setParams()
        self.setSym()

        self.convertUnits()

    def setDS(self):

        with open('ds.csv', 'r') as file:
            dataset = csv.DictReader(file)

            for row in dataset:
                self.atoms.append({'step': int(row['step']),'id': int(row['id']), 'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z']), 'energy': float(row['energy']), 's2': 0.0})

    def setParams(self):

        list = [item['id'] for item in self.atoms]
        self.atmax = max(list)
        list = [item['step'] for item in self.atoms]
        self.stepmax = max(list)

    def setSym(self):
        for i in range(1,self.stepmax):
            lista = list(filter(lambda u: u['step'] == i, self.atoms))
            for atom in lista:
                atom = Symm2D(lista, atom, self.eta_prm, self.rs_prm)
