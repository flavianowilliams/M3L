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

    def __str__(self):
        return (
            f"""Célula unitária ortorrômbica\nVolume: {self.getVolume()['value']*self.getAconv()**3} {self.getVolume()['unit']}\nTotal: {self.getNatoms()} átomos.
            """
        )

class Training(DataSet):

    def __str__(self):
        return(
            f"Total steps: {self.stepmax}\nTotal of atoms: {self.atmax}"
        )