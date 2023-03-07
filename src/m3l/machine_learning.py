import csv
from scipy import constants as cte

class DataSet():

    def __init__(self):

        self.atmax = 0
        self.stepmax = 0
        self.data_sample = list()
        self.setDS()
        self.setParams()

    def setDS(self):

        with open('ds.csv', 'r') as file:
            dataset = csv.DictReader(file)

            for row in dataset:
                self.data_sample.append({'step': int(row['step']), 'id': int(row['id']), 'x': float(row['x']), 'y': float(row['y']), 'z': float(row['z']), 'energy': float(row['energy'])})

    def setParams(self):

        list = [item['id'] for item in self.data_sample]
        self.atmax = max(list)
        list = [item['step'] for item in self.data_sample]
        self.stepmax = max(list)