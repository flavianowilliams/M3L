import csv
from m3l.utils import Conversion

class Atom(Conversion):

    atoms = []

    def __init__(self):

        super().__init__()

    def setMass(self, atm):

        if atm == 'H':
            return 1.008
        elif atm == 'Ar':
            return 39.948

class System(Atom):

    temperature = 0.0
    pressure = 0.0

    def __init__(self, acell, bcell, ccell, filename):

        self.acell = acell
        self.bcell = bcell
        self.ccell = ccell
        self.filename = filename
        self.reload = False

    def setAtoms(self):
        if self.reload:
            with open('system.csv', newline='') as xyz_file:
                file = csv.DictReader(xyz_file)
                for row in file:
                    self.atoms.append({
                        'id': int(row['id']),
                        'symbol': row['symbol'],
                        'mass': float(row['mass']),
                        'x': float(row['x']),
                        'y': float(row['y']),
                        'z': float(row['z']),
                        'vx': float(row['vx']),
                        'vy': float(row['vy']),
                        'vz': float(row['vz']),
                        'fx': float(row['fx']),
                        'fy': float(row['fy']),
                        'fz': float(row['fz']),
                        's2': float(row['s2'])
                        })
        else:
            with open(self.filename, 'r') as xyz_file:
                natoms = xyz_file.readline()
                natoms = int(natoms)
                xyz_file.readline()
                for i in range(natoms):
                    p1, p2, p3, p4 = xyz_file.readline().split(maxsplit=3)
                    p1 = str(p1)
                    p2 = float(p2)
                    p3 = float(p3)
                    p4 = float(p4)
                    mass = self.setMass(p1)
                    self.atoms.append({
                        'id': i,
                        'symbol': p1,
                        'mass': mass,
                        'x': p2,
                        'y': p3,
                        'z': p4,
                        'vx': 0.0,
                        'vy': 0.0,
                        'vz': 0.0,
                        'fx': 0.0,
                        'fy': 0.0,
                        'fz': 0.0,
                        's2': 0.e0
                        })

        return self.atoms

    def setNAtoms(self):

        self.natoms = len(self.atoms)

        return self.natoms

    def setVolume(self):
        self.volume = {'value': self.acell*self.bcell*self.ccell, 'unit': 'A³'}

        return self.volume

    def ccp(self):
        if self.acell > 0.0 and self.bcell > 0.0 and self.ccell > 0.0:
            for at in self.atoms:
                at['x'] = at['x']-self.acell*int(at['x']/self.acell)
                at['y'] = at['y']-self.bcell*int(at['y']/self.bcell)
                at['z'] = at['z']-self.ccell*int(at['z']/self.ccell)
#        else:
#            logging.critical('Constant lattice must be larger than zero!')
#            sys.exit()

