#from scipy import constants as cte
import numpy as np 

class Constants():

    # physical constants
#    KB = cte.value('Boltzmann constant')
#    N0 = cte.value('Avogadro constant')
#    A0 = cte.value('atomic unit of length')
#    HARTREE = cte.value('atomic unit of energy')
#    ELECTRON_MASS = cte.value('atomic unit of mass')
#    PERIOD_BOHR = cte.value('atomic unit of time')

    PERIOD_BOHR = 2.4188843265864e-17
    ELECTRON_MASS = 9.1093837139e-31
    HARTREE = 4.359744722206e-18
    KB = 1.380649e-23
    N0 = 6.02214076e+23
    A0 = 5.29177210544e-11
#    print(PERIOD_BOHR)
#    print(ELECTRON_MASS)
#    print(HARTREE)
#    print(KB)
#    print(N0)
#    print(A0)
    # Conversion from input units to SI 
    ACONV = 1.0e-10                                                # angstrom => metre 
    ECONV = 6.95e-21                                               # kcal/mol => Joule 
    MCONV = 1.0/(1000*N0)                                          # amu => kilogram 
    TIMECONV = 1.0e-12                                             # picosecond => second 
    PCONV = 101325                                                 # atm => Pascal

    # Conversion from input to atomic unit
    ACONV = np.array(ACONV/A0, dtype = np.float64)                 # metre => a0
    ECONV = np.array(ECONV/HARTREE, dtype = np.float64)            # Joule => Hartree
    MCONV = np.array(MCONV/ELECTRON_MASS, dtype = np.float64)      # kilogram => electron mass 
    TIMECONV = np.array(TIMECONV/PERIOD_BOHR, dtype = np.float64)  # second => electron time revolution
    PCONV = np.array(PCONV/(HARTREE/A0**3), dtype = np.float64)    # Pascal => Hartree/a0**3

    # Other conversion
    TEMPCONV = np.array(KB/HARTREE, dtype = np.float64)            # 
    KB = np.array(1.0, dtype = np.float64)                         # J/K => Hartree/K
    CHG = np.array(1.0, dtype = np.float64)                         # Coulomb => a.u.

    # Inverted conversion

    PCONVINV = np.array(1.0e0/PCONV, dtype = np.float64)           # atm^-1

