module structure_module 
  integer :: natom, nfree
  real(8) :: eta
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: atom
contains
!
!-periodic condition contour
!
  subroutine ccp()

    implicit none
    integer :: i 

    do i = 1, natom
      atom(i, 4) = atom(i, 4)-cell(1)*nint(atom(i, 4)/cell(1))
      atom(i, 5) = atom(i, 5)-cell(2)*nint(atom(i, 5)/cell(2))
      atom(i, 6) = atom(i, 6)-cell(3)*nint(atom(i, 6)/cell(3))
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
