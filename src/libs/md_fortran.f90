module Libs 
  integer :: natom
  integer, allocatable, dimension(:) :: nlist
  integer, allocatable, dimension(:, :) :: ilist
  real(8) :: timestep, sigma, tstat, friction
  real(8) :: ekinetic, energy
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:) :: params
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

      if (allocated(fx).eqv..FALSE.) then
        allocate(fx(natom), fy(natom), fz(natom), ea(natom))
      end if

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
    subroutine bonds
      implicit none
      real(8) :: epot, fr, dx, dy, dz, dr 
      
      dx = rx(2)-rx(1)
      dy = ry(2)-ry(1)
      dz = rz(2)-rz(1)
      call mic(dx, dy, dz)
      dr = sqrt(dx**2+dy**2+dz**2)
      epot = 0.5d0*params(1)*(dr-params(2))**2
      fr = -1.0d0*params(1)*(dr-params(2))/dr
      fx(1) = -fr*dx
      fy(1) = -fr*dy
      fz(1) = -fr*dz
      ea(1) = epot
      fx(2) = +fr*dx
      fy(2) = +fr*dy
      fz(2) = +fr*dz
      ea(2) = epot

    end subroutine bonds
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
      real(8) :: depot, epot, fr, dx, dy, dz, dr 
      
      call neighbour_list

      energy = 0.d0
      do i = 1, natom-1
        do j = 1, nlist(i) 
          nj = ilist(i,j)
          dx = rx(nj)-rx(i)
          dy = ry(nj)-ry(i)
          dz = rz(nj)-rz(i)
          call mic(dx, dy, dz)
          dr = sqrt(dx**2+dy**2+dz**2)
          depot = (params(2)/dr)**6
          epot = 4.0d0*params(1)*(depot-1.0d0)*depot
          fr = 24.d0*params(1)*(2.d0*depot-1.d0)*depot/dr**2
          fx(i) = -fr*dx
          fy(i) = -fr*dy
          fz(i) = -fr*dz
          ea(i) = epot
          fx(nj) = +fr*dx
          fy(nj) = +fr*dy
          fz(nj) = +fr*dz
          ea(nj) = epot
          energy = energy+epot
!          if(dr*0.53.le.3.5)print*, dr*0.53, i, nj, ilist(i, j) 
        end do
      end do 

end subroutine vdw
!
!-Neighbours list subroutine function
!
    subroutine neighbour_list
      implicit none
      integer :: i, j, nx
      real(8) :: dx, dy, dz, dr

      if (allocated(nlist).eqv..FALSE.) then
        allocate(nlist(natom), ilist(natom, natom))
      end if 

      do i = 1, natom
        nx = 1
        do j = 1, natom
          if (i.ne.j)then
            dx = rx(j)-rx(i)
            dy = ry(j)-ry(i)
            dz = rz(j)-rz(i)
            call mic(dx, dy, dz)
            dr = sqrt(dx**2+dy**2+dz**2)
            if (dr.le.params(3)) then
              ilist(i, nx) = j
              nx = nx+1
            end if 
          end if 
!          if(i.eq.ilist(i, nx-1))print*, i, ilist(i, nx-1)
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
        rx(i) = rx(i)-cell(1)*int(rx(i)/cell(1))
        ry(i) = ry(i)-cell(2)*int(ry(i)/cell(2))
        rz(i) = rz(i)-cell(3)*int(rz(i)/cell(3))
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

    end subroutine nvt

    subroutine friction_func()
      implicit none

      friction = friction+0.25d0*timestep*(ekinetic-sigma)/(sigma*tstat**2)

    end subroutine friction_func

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

    end subroutine nvt_hoover

end module Libs
