module prepare
    contains
    subroutine teste(filename, at, x, y, z)

        implicit none
    
        integer i
        character(*), intent(in) :: filename
        character(2), intent(out) :: at(2)
        real(8), intent(out) :: x(2), y(2), z(2)

        open(unit=5, file=filename, status='old')

        read(5,*)
        read(5,*)

        do i=1,2
            read(5,*)at(i), x(i), y(i), z(i)
        end do

    end subroutine teste
end module prepare
