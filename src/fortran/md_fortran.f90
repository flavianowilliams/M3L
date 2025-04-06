module Libs 
  integer :: natom
  integer, allocatable, dimension(:) :: nlist
  integer, allocatable, dimension(:,:) :: ilist
  integer, allocatable, dimension(:) :: atp
  real(8) :: timestep, sigma, tstat, pstat, bfactor, press_bath, friction
  real(8) :: ekinetic, energy, temperature, pressure, virial, drmax 
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: params
  real(8), allocatable, dimension(:) :: mass 
  real(8), allocatable, dimension(:) :: rx, ry, rz 
  real(8), allocatable, dimension(:) :: vx, vy, vz 
  real(8), allocatable, dimension(:) :: fx, fy, fz, ea
contains
!
!-Calculate forces a atomic energy
!
  subroutine forces

    implicit none
    integer :: i

    do i = 1, natom
      fx(i) = 0.d0
      fy(i) = 0.d0
      fz(i) = 0.d0
      ea(i) = 0.d0
    end do 

!      call bonds
    call vdw

  end subroutine forces
!
!- intramolecular bonds
!
!  subroutine bonds
!    implicit none
!    real(8) epot, fr, dx, dy, dz, dr 
!    
!    dx = rx(2)-rx(1)
!    dy = ry(2)-ry(1)
!    dz = rz(2)-rz(1)
!    call mic(dx, dy, dz)
!    dr = sqrt(dx**2+dy**2+dz**2)
!    epot = 0.5d0*params(1)*(dr-params(2))**2
!    fr = -1.0d0*params(1)*(dr-params(2))/dr
!    fx(1) = -fr*dx
!    fy(1) = -fr*dy
!    fz(1) = -fr*dz
!    ea(1) = epot
!    fx(2) = +fr*dx
!    fy(2) = +fr*dy
!    fz(2) = +fr*dz
!    ea(2) = epot
!
!  end subroutine bonds
!
!- Bond constraint rules
!
  subroutine mic(dx, dy, dz)

    implicit none
    real(8), intent(inout) :: dx, dy, dz

    dx = dx - cell(1)*nint(dx/cell(1))
    dy = dy - cell(2)*nint(dy/cell(2))
    dz = dz - cell(3)*nint(dz/cell(3))

  end subroutine mic
!
!- Intermolecular interaction
!
  subroutine vdw

    implicit none
    integer :: i, j, nj 
    real(8) :: depot, epot, fr, dx, dy, dz, dr, dvir, encorr, vircorr, drmin, prm1, prm2
    
!    call neighbour_list

    energy = 0.d0
    virial = 0.d0 
    do i = 1, natom-1
      dvir = 0.d0 
      ea(i) = 0.d0
      do j = i+1, natom
        nj = j
        prm1 = sqrt(params(atp(i), 2)*params(atp(nj), 2))
        prm2 = 0.5d0*(params(atp(i), 3)+params(atp(nj), 3))
!      do j = 1, nlist(i) 
!        nj = ilist(i,j)
        dx = rx(nj)-rx(i)
        dy = ry(nj)-ry(i)
        dz = rz(nj)-rz(i)
        call mic(dx, dy, dz)
        dr = sqrt(dx**2.0d0+dy**2.0d0+dz**2.0d0)
        depot = (prm2/dr)**6.0d0 
        epot = 4.0d0*prm1*(depot-1.0d0)*depot
        fr = 24.d0*prm1*(2.d0*depot-1.d0)*depot/dr**2.0d0 
        fx(i) = fx(i)-fr*dx
        fy(i) = fy(i)-fr*dy
        fz(i) = fz(i)-fr*dz
        fx(nj) = fx(nj)+fr*dx
        fy(nj) = fy(nj)+fr*dy
        fz(nj) = fz(nj)+fr*dz
        ea(i) = ea(i)+epot
        dvir = dvir+fr*dr**2.0d0 
        end do
      energy = energy+ea(i)
      virial = virial+dvir
    end do 

    drmin = min(cell(1), cell(2), cell(3))
    drmin = 0.5d0*drmin

!    encorr = 4.0d0*params(1, 3)*(params(1, 4)**12/(9.0d0*drmin**9)-params(1, 4)**6/(3.d0*drmin**3)) 
!    vircorr = -24.0d0*params(1, 3)*(2.0d0*params(1, 4)**12/(9.0d0*drmin**9)-params(1, 4)**6/(3.d0*drmin**3)) 
!
!    encorr = 2.0d0*3.141593d0*natom*natom*encorr/(cell(1)*cell(2)*cell(3))
!    vircorr = 2.0d0*3.141593d0*natom*natom*vircorr/(cell(1)*cell(2)*cell(3))

    encorr = 0.0d0 
    vircorr = 0.0d0 

    energy = energy+encorr
    virial = virial+vircorr

  end subroutine vdw
!
!-Neighbours list subroutine function
!
  subroutine neighbour_list

    implicit none
    integer :: i, j, nx
    real(8) :: dx, dy, dz, dr, drmin

    if (allocated(nlist).eqv..FALSE.) then
      allocate(nlist(natom), ilist(natom, natom))
    end if 

    drmin = min(cell(1), cell(2), cell(3))
    drmin = 0.5*drmin

    do i = 1, natom
      nx = 1
      do j = 1, natom
        if (i.ne.j)then
          dx = rx(j)-rx(i)
          dy = ry(j)-ry(i)
          dz = rz(j)-rz(i)
          call mic(dx, dy, dz)
          dr = sqrt(dx**2+dy**2+dz**2)
          if (dr.lt.drmin) then
            ilist(i, nx) = j
            nx = nx+1
          end if 
        end if 
      end do
      nlist(i) = nx-1
    end do

  end subroutine neighbour_list
!
!-periodic condition contour
!
  subroutine ccp()

    implicit none
    integer :: i 

    do i = 1, natom
      rx(i) = rx(i)-cell(1)*nint(rx(i)/cell(1))
      ry(i) = ry(i)-cell(2)*nint(ry(i)/cell(2))
      rz(i) = rz(i)-cell(3)*nint(rz(i)/cell(3))
    end do 

  end subroutine ccp
!
!-Ensemble nvt berendsen
!
  subroutine ekinetic_func()

    implicit none
    integer :: i 
    real(8) :: sum

    sum = 0.d0
    do i = 1, natom
      sum = sum+mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
    end do

    ekinetic = 0.5d0*sum

  end subroutine ekinetic_func

  subroutine temperature_func
  
    temperature = 2.0d0*ekinetic/(3.0d0*(natom-1))

  end subroutine temperature_func

  subroutine pressure_func

    implicit none 
    real(8) :: volume 
    
    volume = cell(1)*cell(2)*cell(3)

    pressure = (2.0d0*ekinetic+virial)/(3.0d0*volume)

  end subroutine pressure_func
!
!- calculando deslocamento maximo
!
!  subroutine drmax_func
!
!    implicit none
!    integer :: i 
!    real(8) vr 
!
!    drmax = 0.0d0
!    do i = 1, natom
!      vr = sqrt(vx(i)**2+vy(i)**2+vz(i)**2)
!      drmax = max(drmax, vr*timestep)
!    end do
!
!  end subroutine drmax_func
!
!-Ensemble nvt berendsen
!
  subroutine nvt

    implicit none
    integer :: i 
    real(8) :: qui
    
    do i = 1, natom 
      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
      rx(i) = rx(i)+vx(i)*timestep
      ry(i) = ry(i)+vy(i)*timestep
      rz(i) = rz(i)+vz(i)*timestep
    end do 

    call ccp

    call forces

    do i = 1, natom
      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
    end do 

    call ekinetic_func

    qui = sqrt(1.d0+timestep*(sigma/ekinetic-1.0d0)/tstat)

    do i = 1, natom
      vx(i) = vx(i)*qui
      vy(i) = vy(i)*qui
      vz(i) = vz(i)*qui
    end do 

    call ekinetic_func
    call temperature_func

  end subroutine nvt

  subroutine friction_func

    implicit none

    friction = friction+0.25d0*timestep*(ekinetic-sigma)/(sigma*tstat**2)

  end subroutine friction_func

  subroutine eta_func(eta)

    implicit none
    real(8), intent(out) :: eta

    eta = (1.0d0-bfactor*timestep*(press_bath-pressure)/pstat)**(1.0d0/3.0d0)

  end subroutine eta_func

  subroutine nvt_hoover

    implicit none
    integer :: i 
    
    call ekinetic_func
    call friction_func

    do i = 1, natom 
      vx(i) = vx(i)*exp(-0.5d0*timestep*friction)
      vy(i) = vy(i)*exp(-0.5d0*timestep*friction)
      vz(i) = vz(i)*exp(-0.5d0*timestep*friction)
    end do 

    call ekinetic_func
    call friction_func

    do i = 1, natom
      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
      rx(i) = rx(i)+vx(i)*timestep
      ry(i) = ry(i)+vy(i)*timestep
      rz(i) = rz(i)+vz(i)*timestep
    end do 

    call ccp
    call forces

    do i = 1, natom
      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
    end do 

    call ekinetic_func
    call friction_func

    do i = 1, natom
      vx(i) = vx(i)*exp(-0.5d0*timestep*friction)
      vy(i) = vy(i)*exp(-0.5d0*timestep*friction)
      vz(i) = vz(i)*exp(-0.5d0*timestep*friction)
    end do 

    call ekinetic_func
    call friction_func
    call temperature_func

  end subroutine nvt_hoover

  subroutine npt_berendsen

    implicit none
    integer :: i
    real(8) :: eta, qui

    call eta_func(eta)

    do i = 1, natom 
      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
      rx(i) = rx(i)*eta+vx(i)*timestep
      ry(i) = ry(i)*eta+vy(i)*timestep
      rz(i) = rz(i)*eta+vz(i)*timestep
    end do 

    do i = 1,3
      cell(i) = cell(i)*eta 
    end do 

    call ccp
    call forces

    do i = 1, natom 
      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
    end do 

    call ekinetic_func

    qui = sqrt(1.0d0+timestep*(sigma/ekinetic-1.0d0)/tstat)

    do i = 1, natom 
      vx(i) = vx(i)*qui
      vy(i) = vy(i)*qui
      vz(i) = vz(i)*qui
    end do 

    call ekinetic_func
    call temperature_func
    call pressure_func

  end subroutine npt_berendsen

end module Libs
