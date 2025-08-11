import numpy as np
from m3l.structure import System
from m3l.utils import Constants
from m3l import libens, libtherm, libstruct, libforce

class Ensemble(Constants):

    def __init__(self, timestep, force_field, temp_bath = None, press_bath = None, tstat = None, pstat = None, bfactor = None):

        self.timestep = np.array(timestep*self.TIMECONV, dtype = np.float64)
        self.dtimestep = np.array(timestep*self.TIMECONV, dtype = np.float64)
        self.friction = 0.e0

        if tstat:
            self.tstat = np.array(tstat*self.TIMECONV, dtype = np.float64)
        else:
            self.tstat = np.array(2.0*self.TIMECONV, dtype = np.float64)

        if pstat:
            self.pstat = np.array(pstat*self.TIMECONV, dtype = np.float64)
        else:
            self.pstat = np.array(2.0*self.TIMECONV, dtype = np.float64)

        if temp_bath:
            self.temp_bath = np.array(temp_bath*self.TEMPCONV, dtype = np.float64)
        else:
            self.temp_bath = np.array([], dtype = np.float64)

        if press_bath:
            self.press_bath = np.array(press_bath*self.PCONV, dtype = np.float64)
        else:
            self.press_bath = np.array([], dtype = np.float64)

        if bfactor:
            self.bfactor = np.array(bfactor, dtype = np.float64)
        else:
            self.bfactor = self.BETAFACTOR
    
        self.bfactor = self.bfactor*self.PCONVINV

        self.setForceField(force_field)

    def setForceField(self, force_field):

        self.nvdw = np.array(force_field[0], dtype=np.int32)
        self.rvdw = np.array(force_field[1], dtype=np.float64)
        self.rcoul = np.array(force_field[2], dtype=np.float64)
        self.ijinter = np.array(force_field[3][0], dtype=np.float64)
        self.prmvdw = np.array(force_field[3][1], dtype=np.float64)

    def nve_verlet(self):

        print("não valendo...")

    def nvt_verlet(self):

        libens.temp_friction = self.sys.temp_friction
        libens.press_friction = self.sys.press_friction

#-atualizando coordenadas canônicas

        libens.tstat = self.tstat
        libens.temp_bath = self.temp_bath
        libens.nfree = self.sys.nfree
        libens.timestep = self.timestep
        libens.natom = self.sys.natom
        libens.ra = self.sys.ra
        libens.va = self.sys.va
        libens.fa = self.sys.fa
        libens.mass = self.sys.mass

        libens.nvt_berendsen_velocity()
        libens.nvt_berendsen_position()

#-aplicando condicoes de contorno periodicas
        
        libstruct.natom = self.sys.natom
        libstruct.cell = self.sys.cell
        libstruct.ra = libens.ra
        libstruct.ccp()

#-calculando campo de força

        libforce.natom = libstruct.natom
        libforce.atype = self.sys.atype
        libforce.charge = self.sys.charge
        libforce.ra = libstruct.ra 
        libforce.fa = libens.fa 
        libforce.ea = self.sys.ea 
        libforce.sites = self.sys.sites
        libforce.nsites = self.sys.nsites
        libforce.ma = self.sys.ma
        libforce.cell = libstruct.cell
        libforce.nvdw = self.nvdw 
        libforce.rvdw = self.rvdw
        libforce.rcoul = self.rcoul
        libforce.ijinter = self.ijinter
        libforce.prmvdw = self.prmvdw
        libforce.setforces()

#-atualizando coordenadas canônicas

        libens.fa = libforce.fa
        libens.nvt_berendsen_velocity()

#-calculando energia cinética

        libtherm.natom = libens.natom
        libtherm.va = libens.va 
        libtherm.mass = libens.mass
        libtherm.setekinetic()

#-atualizando coordenadas canônicas

        libens.ekinetic = libtherm.ekinetic
        libens.nvt_berendsen_velocity_update()

#-calculando energia cinética e temperatura

        libtherm.va = libens.va 
        libtherm.mass = libens.mass
        libtherm.nfree = libens.nfree
        libtherm.setekinetic()
        libtherm.settemperature()

#-calculando pressão 

        libtherm.volume = self.sys.volume
        libtherm.virial = libforce.virial
        libtherm.setpressure()

    def nvt_hoover(self):

#- calculando parametro xi

        libens.ekinetic = self.sys.ekinetic
        libens.tstat = self.tstat
        libens.temp_bath = self.temp_bath
        libens.nfree = self.sys.nfree
        libens.timestep = self.timestep
        libens.natom = self.sys.natom
        libens.ra = self.sys.ra
        libens.va = self.sys.va
        libens.fa = self.sys.fa
        libens.mass = self.sys.mass
        libens.temp_friction = self.sys.temp_friction

        libens.settempfrictionhoover()

#-atualizando coordenadas canônicas

        libens.nvt_hoover_velocity()

#-calculando energia cinética

        libtherm.natom = self.sys.natom
        libtherm.va = libens.va 
        libtherm.mass = libens.mass
        libtherm.setekinetic()

#-atualizando coordenadas canônicas

        libens.ekinetic = libtherm.ekinetic
        libens.settempfrictionhoover()
        libens.nvt_hoover_velocity_update()
        libens.nvt_hoover_position()

#-aplicando condicoes de contorno periodicas
        
        libstruct.natom = self.sys.natom
        libstruct.cell = self.sys.cell
        libstruct.ra = libens.ra
        libstruct.ccp()

#-calculando campo de força

        libforce.natom = libstruct.natom
        libforce.atype = self.sys.atype
        libforce.charge = self.sys.charge
        libforce.ra = libstruct.ra 
        libforce.fa = self.sys.fa 
        libforce.ea = self.sys.ea 
        libforce.sites = self.sys.sites
        libforce.nsites = self.sys.nsites
        libforce.ma = self.sys.ma
        libforce.cell = libstruct.cell
        libforce.nvdw = self.nvdw 
        libforce.rvdw = self.rvdw
        libforce.rcoul = self.rcoul
        libforce.ijinter = self.ijinter
        libforce.prmvdw = self.prmvdw
        libforce.setforces()

#-atualizando coordenadas canônicas

        libens.fa = libforce.fa
        libens.nvt_hoover_velocity_update()

#-calculando energia cinética

        libtherm.va = libens.va 
        libtherm.setekinetic()

#-atualizando coordenadas canônicas

        libens.ekinetic = libtherm.ekinetic
        libens.settempfrictionhoover()

        libens.nvt_hoover_velocity()

#-calculando energia cinética

        libtherm.va = libens.va 
        libtherm.setekinetic()

#-atualizando coeficiente de friccao

        libens.ekinetic = libtherm.ekinetic
        libens.settempfrictionhoover()
        
#-calculando energia cinética e temperatura

        libtherm.va = libens.va 
        libtherm.nfree = self.sys.nfree 
        libtherm.setekinetic()
        libtherm.settemperature()

#-calculando pressão 

        libtherm.virial = libforce.virial
        libtherm.cell = self.sys.cell
        libtherm.setvolume()
        libtherm.setpressure()

    def npt_hoover(self):

#-atribuindo valores fixos

        libens.tstat = self.tstat
        libens.pstat = self.pstat
        libens.temp_bath = self.temp_bath
        libens.press_bath = self.press_bath
        libens.nfree = self.sys.nfree
        libens.timestep = self.timestep
        libens.natom = self.sys.natom
        libtherm.natom = self.sys.natom
        libforce.natom = self.sys.natom
        libens.mass = self.sys.mass
        libstruct.mass = self.sys.mass
        libtherm.mass = self.sys.mass
        libens.temp_friction = self.sys.temp_friction
        libens.press_friction = self.sys.press_friction

#-calculando energia cinética

        libtherm.va = self.sys.va 
        libtherm.setekinetic()

#-atualizando coeficiente de friccao

        self.setTempFriction()

#-atualizando pressão e volume

        libtherm.pressure = self.sys.pressure
        libtherm.volume = self.sys.volume

#-atualizando coeficiente de friccao

        self.setPressFriction()

#-atualizando coordenadas canônicas

        lx = np.exp(-0.5e0*self.timestep*libens.press_friction)
        for i in range(self.sys.natom):
            self.sys.va[i][0] = self.sys.va[i][0]*lx
            self.sys.va[i][1] = self.sys.va[i][1]*lx
            self.sys.va[i][2] = self.sys.va[i][2]*lx

#-calculando energia cinética

        libtherm.va = self.sys.va
        libtherm.setekinetic()

#-calculando pressão 

        libtherm.virial = self.sys.virial
        libtherm.setpressure()

#-atualizando coeficiente de friccao

        self.setPressFriction()

#-atualizando coeficiente de friccao

        self.setTempFriction()

#-atualizando coordenadas canonicas

        for i in range(self.sys.natom):
            for j in range(3):
                self.sys.va[i][j] = self.sys.va[i][j]+0.5e0*self.timestep*self.sys.fa[i][j]/self.sys.mass[i]

#-atualizando parâmetros de rede e volume da supercélula

        libstruct.eta = np.exp(libens.press_friction*self.timestep)
        libstruct.cell = self.sys.cell
        libstruct.setcell()

        libtherm.cell = libstruct.cell
        libtherm.setvolume()

#-aplicando condicoes de contorno periodicas
        
        libstruct.natom = self.sys.natom
        libstruct.ra = self.sys.ra
        libstruct.ccp()

#-calculando centro de massa do sistema

        libstruct.setrcm()

#-atualizando coordenadas canonicas

        lx = libstruct.eta
        for i in range(self.sys.natom):
            for j in range(3):
                self.sys.ra[i][j] = (self.sys.ra[i][j]-libstruct.rcm[j])*lx+self.timestep*self.sys.va[i][j]+libstruct.rcm[j]

#-calculando campo de força

        self.setForces()

#-atualizando coordenadas canonicas

        for i in range(self.sys.natom):
            for j in range(3):
                self.sys.va[i][j] = self.sys.va[i][j]+0.5e0*self.timestep*self.sys.fa[i][j]/self.sys.mass[i]

#-atualizando coeficiente de friccao

        self.setTempFriction()

#-atualizando coeficiente de friccao

        self.setPressFriction()

#-atualizando coordenadas canônicas

        lx = np.exp(-0.5e0*self.timestep*libens.press_friction)
        for i in range(self.sys.natom):
            self.sys.va[i][0] = self.sys.va[i][0]*lx
            self.sys.va[i][1] = self.sys.va[i][1]*lx
            self.sys.va[i][2] = self.sys.va[i][2]*lx

#-calculando energia cinética

        libtherm.va = self.sys.va
        libtherm.setekinetic()

#-calculando pressão 

        libtherm.virial = libforce.virial
        libtherm.setpressure()

#-atualizando coeficiente de friccao

        self.setPressFriction()

#-atualizando coeficiente de friccao

        self.setTempFriction()

#-atualizando coordenadas canônicas

        lx = np.exp(-0.5e0*self.timestep*libens.press_friction)
        for i in range(self.sys.natom):
            self.sys.va[i][0] = self.sys.va[i][0]*lx
            self.sys.va[i][1] = self.sys.va[i][1]*lx
            self.sys.va[i][2] = self.sys.va[i][2]*lx

#-atualizando coeficiente de friccao

        self.setPressFriction()

#-atualizando coeficiente de friccao

        self.setTempFriction()

#-calculando momento do centro de massa do sistema

        libstruct.va = self.sys.va
        libstruct.setvcm()

#-atualizando coordenadas canônicas

        vcm = libstruct.vcm
        for i in range(self.sys.natom):
            self.sys.va[i][0] = self.sys.va[i][0]-vcm[0]
            self.sys.va[i][1] = self.sys.va[i][1]-vcm[1]
            self.sys.va[i][2] = self.sys.va[i][2]-vcm[2]

#-calculando energia cinética e temperatura

        libtherm.va = self.sys.va
        libtherm.setekinetic()
        libtherm.nfree = self.sys.nfree
        libtherm.settemperature()

#-calculando pressão 

        libtherm.virial = libforce.virial
        libtherm.setpressure()

#-atribuindo variaveis canonicas para ser utilizado posteriormente

        libens.ra = self.sys.ra 
        libens.va = self.sys.va 

    def setTempFriction(self):

#-atualizando coeficiente de friccao

        libens.ekinetic = libtherm.ekinetic
        libens.settempfrictionhoovernpt()

#-atualizando coordenadas canônicas

#        libens.npt_hoover_velocity()
        lx = np.exp(-0.25e0*self.timestep*libens.temp_friction)
        for i in range(self.sys.natom):
            self.sys.va[i][0] = self.sys.va[i][0]*lx
            self.sys.va[i][1] = self.sys.va[i][1]*lx
            self.sys.va[i][2] = self.sys.va[i][2]*lx

#-calculando energia cinética

        libtherm.va = self.sys.va 
        libtherm.setekinetic()

#-atualizando coeficiente de friccao

        libens.ekinetic = libtherm.ekinetic
        libens.settempfrictionhoovernpt()

    def setPressFriction(self):

#- calculando coeficiente de friccao
        
        libens.pressure = libtherm.pressure
        libens.volume = libtherm.volume
        libens.setpressfrictionhoover()

    def setForces(self):

        libforce.atype = self.sys.atype
        libforce.charge = self.sys.charge
        libforce.ra = libstruct.ra 
        libforce.fa = self.sys.fa 
        libforce.ea = self.sys.ea 
        libforce.sites = self.sys.sites
        libforce.nsites = self.sys.nsites
        libforce.ma = self.sys.ma
        libforce.cell = libstruct.cell
        libforce.nvdw = self.nvdw 
        libforce.rvdw = self.rvdw
        libforce.rcoul = self.rcoul
        libforce.ijinter = self.ijinter
        libforce.prmvdw = self.prmvdw
        libforce.setforces()

#    def npt_verlet(self):
#
#        libens.temp_friction = self.sys.temp_friction
#        libens.press_friction = self.sys.press_friction
#
##- calculando parametro eta
#        #
#        libens.bfactor = self.bfactor
#        libens.timestep = self.timestep
#        libens.press_bath = self.press_bath
#        libens.pressure = self.sys.pressure
#        libens.pstat = self.pstat
#
#        libens.seteta()
#
##-atualizando coordenadas canônicas
#
#        libens.natom = self.sys.natom
#        libens.ra = self.sys.ra 
#        libens.va = self.sys.va
#        libens.fa = self.sys.fa 
#        libens.mass = self.sys.mass
#        libens.npt_berendsen_velocity()
#        libens.npt_berendsen_position()
#
##-atualizando parâmetros de rede e volume da supercélula
#
#        libstruct.eta = libens.eta
#        libstruct.cell = self.sys.cell
#        libstruct.setcell()
#
#        libtherm.cell = libstruct.cell
#        libtherm.setvolume()
#
##-atualizando coordenadas atômicas
#        
#        libstruct.natom = self.sys.natom
#        libstruct.ra = libens.ra
#        libstruct.ccp()
#
##-calculando campo de força
#
#        libforce.natom = libstruct.natom
#        libforce.atype = self.sys.atype
#        libforce.charge = self.sys.charge
#        libforce.ra = libstruct.ra 
#        libforce.fa = self.sys.fa 
#        libforce.ea = self.sys.ea 
#        libforce.sites = self.sys.sites
#        libforce.nsites = self.sys.nsites
#        libforce.ma = self.sys.ma
#        libforce.cell = libstruct.cell
#        libforce.nvdw = self.nvdw 
#        libforce.rvdw = self.rvdw
#        libforce.rcoul = self.rcoul
#        libforce.ijinter = self.ijinter
#        libforce.prmvdw = self.prmvdw
#        libforce.setforces()
#
##-atualizando coordenadas canônicas
#
#        libens.fa = libforce.fa
#        libens.npt_berendsen_velocity()
#        
##-calculando energia cinética
#
#        libtherm.natom = self.sys.natom
#        libtherm.va = libens.va 
#        libtherm.mass = libens.mass
#        libtherm.setekinetic()
#
##-atualizando coordenadas canônicas
#
#        libens.nfree = self.sys.nfree
#        libens.ekinetic = libtherm.ekinetic
#        libens.tstat = self.tstat
#        libens.temp_bath = self.temp_bath
#        libens.npt_berendsen_velocity_update()
#        
##-calculando energia cinética e temperatura
#
#        libtherm.va = libens.va
#        libtherm.setekinetic()
#        libtherm.nfree = self.sys.nfree
#        libtherm.settemperature()
#
##-calculando pressão 
#
#        libtherm.virial = libforce.virial
#        libtherm.setpressure()
#
    def hook(self, system):

        self.sys = system

    def hook_output(self):

        description = self.sys.description

        self.sys = System()

        self.sys.description = description
        self.sys.cell = libstruct.cell
        self.sys.volume = libtherm.volume
        self.sys.natom = libens.natom
        self.sys.nfree = libens.nfree
        self.sys.sites = libforce.sites
        self.sys.nsites = libforce.nsites
        self.sys.atype = libforce.atype
        self.sys.charge = libforce.charge
        self.sys.mass = libens.mass
        self.sys.ra = libens.ra
        self.sys.va = libens.va
        self.sys.ma = libforce.ma
        self.sys.fa = libforce.fa 
        self.sys.ea = libforce.ea 
#        
        self.sys.temperature = libtherm.temperature
        self.sys.press_friction = libens.press_friction
        self.sys.temp_friction = libens.temp_friction
        self.sys.pressure = libtherm.pressure
        self.sys.ekinetic = libtherm.ekinetic
        self.sys.epotential = libforce.energy
        self.sys.virial = libforce.virial

        return self.sys

    def __call__(self, system):
        self.hook(system)
        if self.temp_bath.size > 0:
            if self.press_bath.size > 0:
                self.npt_hoover()
            else:
                self.nvt_hoover()
        else:
            self.nve_verlet()
        output = self.hook_output()
        return output 

