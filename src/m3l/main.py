from m3l.machine_learning import DataSet
from m3l.molecular_dynamics import System

class Structure(System):

    def __init__(self, a, b, c, filename, eta_prm, rs_prm):

        self.setAcell(a)
        self.setBcell(b)
        self.setCcell(c)
        self.eta_prm = eta_prm
        self.rs_prm = rs_prm
        self.setXYZ(filename)
        self.setCCP()

        self.convertUnits()

        self.setVolume()

    def convertUnits(self):

        self.setAcell(self.getAcell()/self.getAconv())
        self.setBcell(self.getBcell()/self.getAconv())
        self.setCcell(self.getCcell()/self.getAconv())
        self.eta_prm=self.eta_prm/(1/self.getAconv()**2)
        self.rs_prm=self.rs_prm/self.getAconv()

        for at in self.atoms:
            at['x'] = at['x']/self.getAconv()
            at['y'] = at['y']/self.getAconv()
            at['z'] = at['z']/self.getAconv()

    def __str__(self):

        return (
            f"""Célula unitária ortorrômbica
            \nVolume: {self.getVolume()['value']*self.getAconv()**3}
            {self.getVolume()['unit']}
            \nTotal: {self.getNatoms()} átomos.
            """
        )

class Training(DataSet):

    def __init__(self, eta_prm, rs_prm):

        self.eta_prm = eta_prm
        self.rs_prm = rs_prm
        self.atoms = list()
        self.atmax = 0
        self.stepmax = 0
        self.setDS()
        self.setParams()
        self.setSym()

        self.convertUnits()

    def convertUnits(self):

        self.eta_prm=self.eta_prm/(1/self.getAconv()**2)
        self.rs_prm=self.rs_prm/self.getAconv()

        for at in self.atoms:
            at['x'] = at['x']/self.getAconv()
            at['y'] = at['y']/self.getAconv()
            at['z'] = at['z']/self.getAconv()

    def __str__(self):

        return (f"Total steps: {self.stepmax}\nTotal of atoms: {self.atmax}")
