from scipy import constants as cte
import numpy as np 

class Constants():

    # physical constants
    KB = cte.value('Boltzmann constant')
    N0 = cte.value('Avogadro constant')
    A0 = cte.value('atomic unit of length')
    HARTREE = cte.value('atomic unit of energy')
    ELECTRON_MASS = cte.value('atomic unit of mass')
    PERIOD_BOHR = cte.value('atomic unit of time')

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
    TEMPCONV = np.array(KB/HARTREE, dtype = np.float64)            # J/K => Hartree/K
    KB = np.array(1.0, dtype = np.float64)                         # J/K => Hartree/K

    # Inverted conversion

    PCONVINV = np.array(1.0e0/PCONV, dtype = np.float64)           # atm^-1

    # Some important constants

    BETAFACTOR = np.array(5.0e-5, dtype = np.float64)              # Isothermal compressibility of liquid water (atm)
