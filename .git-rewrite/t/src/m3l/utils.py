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

    # Conversion from SI to input units
    ACONV = 1.0e10                             # metre => angstrom
    ECONV = 1.438689e20                        # Joule => kcal/mol
    MCONV = 1000*N0                            # kilogram => molar mass
    TIMECONV = 1.0e12                          # second => picosecond

    # Conversion from input to atomic unit
    ACONV = np.array(ACONV*A0)                 # angstrom => a0
    ECONV = np.array(ECONV*HARTREE)            # kcal/mol => Hartree
    MCONV = np.array(MCONV*ELECTRON_MASS)      # molar mass => amu
    TIMECONV = np.array(TIMECONV*PERIOD_BOHR)  # picosecond => electron time revolution

    # Other conversion
    KB = np.array(KB/HARTREE)                  # J/K => Hartree/K
    TEMPCONV = np.array(1.0e0)
#    TEMPCONV = np.array(HARTREE/KB)

