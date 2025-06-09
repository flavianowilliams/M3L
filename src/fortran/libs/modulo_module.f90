module fortran_module 
  integer :: natom, nfree
  real(8), allocatable, dimension(:, :) :: atom
contains
!
!- preparing variables
!
  subroutine teste()

    print*, natom

  end subroutine teste

end module fortran_module
