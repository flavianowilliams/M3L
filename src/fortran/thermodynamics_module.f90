module thermodynamics_module 
  integer :: natom, nfree
  real(8) :: ekinetic, temperature, pressure, virial, volume, eta 
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: atom
contains
!
!-Update volume 
!
  subroutine setVolume()

    volume = cell(1)*cell(2)*cell(3)    
    
  end subroutine setVolume
!
  subroutine setVolumex()

    volume = volume*eta**3.0d0
    
  end subroutine setVolumex
!
!-Ensemble nvt berendsen
!
  subroutine setEkinetic()

    implicit none
    integer :: i 
    real(8) :: sum

    sum = 0.d0
    do i = 1, natom
      sum = sum+atom(i, 2)*(atom(i, 7)**2+atom(i, 8)**2+atom(i, 9)**2)
    end do

    ekinetic = 0.5d0*sum

  end subroutine setEkinetic

  subroutine setTemperature
  
    temperature = 2.0d0*ekinetic/nfree

  end subroutine setTemperature

  subroutine setPressure

    implicit none 
    
    pressure = (2.0d0*ekinetic+virial)/(3.0d0*volume)

  end subroutine setPressure

end module thermodynamics_module
