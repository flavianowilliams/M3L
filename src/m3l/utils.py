from scipy import constants as cte

class Constant():
    __aconv = 1.e+10*cte.value('atomic unit of length')

    def getAconv(self):
        return self.__aconv

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
