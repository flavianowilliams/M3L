from m3l.structure import System
from m3l.symmetry import Symm2D

class Structure(Symm2D,System):

    def __init__(self, a, b, c, filename, prm, rs):

        self.setAcell(a)
        self.setBcell(b)
        self.setCcell(c)
        self.setXYZ(filename)
        self.eta_prm = prm
        self.rs_prm = rs

        System.__init__(self, a, b, c, filename)
        Symm2D.__init__(self, self.atoms, prm, rs)

    def __str__(self):
        return (
            '\nCélula unitária ortorrômbica\n\n'
            +'Volume: {} {}\n\n'.format(self.getVolume()['value']*self.getAconv()**3, self.getVolume()['unit'])
            +'Total: {} átomos\n'.format(self.getNatoms())
        )
