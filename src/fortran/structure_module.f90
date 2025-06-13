module structure_module 
  integer :: natom
  real(8) :: eta
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: ra 
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

end module structure_module
