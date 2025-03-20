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
    ACONV = 1.0e-10                            # angstrom => metre 
    ECONV = 6.95e-21                           # kcal/mol => Joule 
    MCONV = 1.0/(1000*N0)                      # amu => kilogram 
    TIMECONV = 1.0e-12                         # picosecond => second 
    PCONV = 101325                             # atm => Pascal

    # Conversion from input to atomic unit
    ACONV = np.array(ACONV/A0)                 # metre => a0
    ECONV = np.array(ECONV/HARTREE)            # Joule => Hartree
    MCONV = np.array(MCONV/ELECTRON_MASS)      # kilogram => electron mass 
    TIMECONV = np.array(TIMECONV/PERIOD_BOHR)  # second => electron time revolution
    PCONV = np.array(PCONV/(HARTREE/A0**3))    # Pascal => Hartree/a0**3

    # Other conversion
    KB = np.array(KB/HARTREE)                  # J/K => Hartree/K
    TEMPCONV = KB 
    KB = np.array(1.0)                  # J/K => Hartree/K

    # Inverted conversion

    PCONVINV = 1.0e0/PCONV                     # atm^-1

    # Some important constants

    BETAFACTOR = 5.0e-5                        # Isothermal compressibility of liquid water (atm)
