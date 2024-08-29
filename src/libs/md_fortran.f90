module libs 
  integer :: natom
  integer, allocatable, dimension(:) :: nlist
  integer, allocatable, dimension(:, :) :: ilist
  real(8) :: timestep
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:) :: params
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
      integer :: bnd
      real(8) :: epot, fr, dx, dy, dz, dr 
      
      dx = rx(2)-rx(1)
      dy = ry(2)-ry(1)
      dz = rz(2)-rz(1)
      call bond_constraint(dx, dy, dz)
      dr = sqrt(dx**2+dy**2+dz**2)
      epot = 0.5d0*params(1)*(dr-params(2))**2
      fr = -1.0d0*params(1)*(dr-params(2))
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
    subroutine bond_constraint(dx, dy, dz)
      implicit none
      real(8), intent(inout) :: dx, dy, dz

      dx = dx - cell(1)*int(2.0d0*dx/cell(1))
      dy = dy - cell(2)*int(2.0d0*dy/cell(2))
      dz = dz - cell(3)*int(2.0d0*dz/cell(3))

    end subroutine bond_constraint
!
!- Intermolecular interaction
!
    subroutine vdw
      implicit none
      integer :: i, j 
      real(8) :: depot, epot, fr, dx, dy, dz, dr, vel_max
      
!      vel_max = 0.d0
!      do i = 1, natom 
!        vel_max = max(vel_max, sqrt(vx(i)**2+vy(i)**2+vz(i)**2)
!      end do 
!      if (vel_max*timestep.ge.params(3)) then
        call neighbour_list
!      end if 

      do i = 1, natom-1
        do j = 1, nlist(i)
          dx = rx(ilist(i,j))-rx(i)
          dy = ry(ilist(i,j))-ry(i)
          dz = rz(ilist(i,j))-rz(i)
          dr = sqrt(dx**2+dy**2+dz**2)
          depot = (params(2)/dr)**6
          epot = 4.0d0*params(1)*(depot-1.0d0)*depot
          fr = 24.d0*params(1)*(2.d0*depot-1.d0)*depot/dr**2
          fx(i) = -fr*dx
          fy(i) = -fr*dy
          fz(i) = -fr*dz
          ea(i) = epot
          fx(j) = +fr*dx
          fy(j) = +fr*dy
          fz(j) = +fr*dz
          ea(j) = epot
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
        allocate(nlist(natom-1), ilist(natom-1, natom-1))
      end if 

      do i = 1, natom-1
        nx = 1
        do j = i+1, natom
          dx = rx(j)-rx(i)
          dy = ry(j)-ry(i)
          dz = rz(j)-rz(i)
          dr = sqrt(dx**2+dy**2+dz**2)
          if (dr.le.params(3)) then
            ilist(i, nx) = j
            nx = nx+1
          end if 
        end do
        nlist(i) = nx-1
      end do

    end subroutine neighbour_list
end module libs
