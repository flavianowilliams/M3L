module forcefield
    contains
    subroutine neighbour_list(rx, ry, rz, rcutoff, natom, ilist)

        implicit none
    

        integer :: i,j,nx
        real(8) :: dr
        integer, intent(in) :: natom
        integer, intent(out) :: ilist(natom, natom)
        real(8), intent(in) :: rx(natom), ry(natom), rz(natom), rcutoff

        do i=1,natom
            nx=1
            do j=i+1,natom
                dr = sqrt((rx(j)-rx(i))**2+(ry(j)-ry(i))**2+(rz(j)-rz(i))**2)
                if (dr.le.rcutoff)then
                    ilist(i, nx) = j
                    nx=nx+1
                end if
            end do
        end do

    end subroutine neighbour_list
end module forcefield
