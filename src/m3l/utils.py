from scipy import constants as cte

class Constant():
    __aconv = 1.e+10*cte.value('atomic unit of length')

    def getAconv(self):
        return self.__aconv

