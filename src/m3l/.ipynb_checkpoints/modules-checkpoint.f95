module module
    implicit none
    contains
    subroutine readxyz(input, output)

        character(*), intent(in) :: input
        integer, intent(out) :: output

!    open(unit=1, file=filename, status="old", action="read")

        output = 100

    end subroutine readxyz
end module module
