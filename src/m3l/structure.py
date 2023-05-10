from m3l.utils import Conversion

class Atom(Conversion):

    atoms = []

    def __init__(self):

        super().__init__()

class System(Atom):

    temperature = 0.0
    pressure = 0.0

    def __init__(self, acell, bcell, ccell, filename):

        self.acell = acell
        self.bcell = bcell
        self.ccell = ccell
        self.filename = filename

    def setAtoms(self):
        with open(self.filename, 'r') as xyz_file:
            natoms = xyz_file.readline()
            natoms = int(natoms)
            xyz_file.readline()
            for i in range(natoms):
                p1, p2, p3, p4 = xyz_file.readline().split(maxsplit=3)
                p2 = float(p2)
                p3 = float(p3)
                p4 = float(p4)
                self.atoms.append({'id': i, 'atom': p1, 'x': p2, 'y': p3, 'z': p4, 'energy': 0.e0, 's2': 0.e0})

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

