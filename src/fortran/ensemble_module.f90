module ensemble_module 
  integer :: natom, nfree
  real(8) :: timestep, tstat, pstat, bfactor, temp_bath, press_bath, friction, eta
  real(8) :: ekinetic, energy, temperature, pressure, virial, temp_friction, press_friction, volume 
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: atom
contains
!
!-Ensemble nvt berendsen
!
  subroutine setEta()

    implicit none

    eta = (1.0d0-bfactor*timestep*(press_bath-pressure)/pstat)**(1.0d0/3.0d0)

  end subroutine setEta

  subroutine temperature_func()
  
    temperature = 2.0d0*ekinetic/nfree

  end subroutine temperature_func

  subroutine pressure_func()

    implicit none 
    
    pressure = (2.0d0*ekinetic+virial)/(3.0d0*volume)

  end subroutine pressure_func

  subroutine friction_func()

    implicit none

    real(8) :: sigma

    sigma = 0.5d0*nfree*temp_bath
    friction = friction+0.25d0*timestep*(ekinetic-sigma)/(sigma*tstat**2)

  end subroutine friction_func

  subroutine setTempFriction()

    implicit none

    real(8) :: term_mass, bar_mass, func

    term_mass = nfree*temp_bath*tstat**2.d0 
    bar_mass = nfree*temp_bath*pstat**2.0d0

    func = 0.5d0*timestep*(nfree*(temperature-temp_bath)+(bar_mass*press_friction**2.0d0-temp_bath))/term_mass

    temp_friction = temp_friction+func

  end subroutine setTempFriction

  subroutine press_friction_func

    implicit none

    real(8) :: term_mass, bar_mass, func

    term_mass = nfree*temp_bath*tstat**2.d0 
    bar_mass = nfree*temp_bath*pstat**2.0d0

    func = 0.5d0*timestep*(3.0d0*volume*(pressure-press_bath)/bar_mass-temp_friction*press_friction)

    press_friction = press_friction + func

  end subroutine press_friction_func

  subroutine npt_berendsen_a()

    implicit none

    integer :: i

    do i = 1, natom 
      atom(i, 7) = atom(i, 7)+atom(i, 10)*0.5d0*timestep/atom(i, 2)
      atom(i, 8) = atom(i, 8)+atom(i, 11)*0.5d0*timestep/atom(i, 2)
      atom(i, 9) = atom(i, 9)+atom(i, 12)*0.5d0*timestep/atom(i, 2)
    end do 

  end subroutine npt_berendsen_a

  subroutine npt_berendsen_b()

    implicit none

    integer :: i

    do i = 1, natom 
      atom(i, 4) = atom(i, 4)*eta+atom(i, 7)*timestep
      atom(i, 5) = atom(i, 5)*eta+atom(i, 8)*timestep
      atom(i, 6) = atom(i, 6)*eta+atom(i, 9)*timestep
    end do 

  end subroutine npt_berendsen_b

  subroutine npt_berendsen_c()

    implicit none

    integer :: i
    real(8) :: sigma, qui

    sigma = 0.5d0*nfree*temp_bath
    qui = sqrt(1.0d0+timestep*(sigma/ekinetic-1.0d0)/tstat)

    do i = 1, natom 
      atom(i, 7) = atom(i, 7)*qui
      atom(i, 8) = atom(i, 8)*qui
      atom(i, 9) = atom(i, 9)*qui
    end do 

  end subroutine npt_berendsen_c

end module ensemble_module
