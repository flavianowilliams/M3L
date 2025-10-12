module structure_module 
  integer :: natom
  real(8) :: eta
  real(8), dimension(3) :: cell, rcm, vcm
  real(8), allocatable, dimension(:) :: mass
  real(8), allocatable, dimension(:, :) :: ra, va 
  real(8), allocatable, dimension(:, :) :: atom
contains
!
!-periodic condition contour
!
  subroutine ccp()

    implicit none
    integer :: i, j 

    do i = 1, natom
      do j = 1, 3
        ra(i, j) = ra(i, j)-cell(j)*nint(ra(i, j)/cell(j))
!        ra(i, j) = ra(i, j)-cell(j)*min(2.0d0, nint(ra(i, j)/cell(j)))
      end do 
    end do 

  end subroutine ccp

  subroutine setCell()

    implicit none
    integer :: i 

    do i = 1,3
      cell(i) = cell(i)*eta 
    end do 

  end subroutine setCell

  subroutine setRcm()

    implicit none

    integer :: i, j
    real(8) :: mtotal

    do i = 1, 3
      rcm(i) = 0.d0
    end do 

    mtotal = 0.d0
    do i = 1, natom
      do j = 1, 3
        rcm(j) = rcm(j)+ra(i, j)*mass(i)
      end do 
      mtotal = mtotal+mass(i)
    end do 

    do i = 1, 3
      rcm(i) = rcm(i)/mtotal
    end do 

  end subroutine setRcm

  subroutine setVcm()

    implicit none

    integer :: i, j
    real(8) :: mtotal

    do i = 1, 3
      vcm(i) = 0.d0
    end do 

    mtotal = 0.d0
    do i = 1, natom
      do j = 1, 3
        vcm(j) = vcm(j)+va(i, j)*mass(i)
      end do 
      mtotal = mtotal+mass(i)
    end do 

    do i = 1, 3
      vcm(i) = vcm(i)/mtotal
    end do 

  end subroutine setVcm

end module structure_module
