from m3l.structure import System
from m3l.symmetry import Symm2D

class Structure(System):

#    def __init__(self, a, b, c, filename):
#
#        self.setAcell(a)
#        self.setBcell(b)
#        self.setCcell(c)
#        self.setXYZ(filename)
#        self.eta_prm = prm
#        self.rs_prm = rs

#        System(self).__init__(a, b, c, filename)

#    self.convertUnits()

#        Symm2D(Structure, self).__init__(self.atoms, prm, rs)
#   
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
            '\nCélula unitária ortorrômbica\n\n'
            +'Volume: {} {}\n\n'.format(self.getVolume()['value']*self.getAconv()**3, self.getVolume()['unit'])
            +'Total: {} átomos\n'.format(self.getNatoms())
        )
