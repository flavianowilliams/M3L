module force_field_module 
  integer :: natom, nvdw, nlist
  integer, allocatable, dimension(:) :: ilist, ma 
  integer, allocatable, dimension(:,:) :: sites 
  real(8) :: virial, rvdw, rcoul, energy 
  integer, allocatable, dimension(:) :: atype, nsites
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: prmvdw, ijinter 
  real(8), allocatable, dimension(:, :) :: ra, fa
  real(8), allocatable, dimension(:) :: ea, charge
contains
!
!-Calculate forces a atomic energy
!
  subroutine setForces()

    implicit none
    integer :: i, j, nx
    real(8) :: drmin

    if (allocated(prmvdw).eqv..FALSE.) then
      allocate(prmvdw(nvdw, 2))
    end if 

    if (allocated(ilist).eqv..FALSE.) then
      nx = size(nsites)
      allocate(ilist(nx))
    end if 

    drmin = 0.5d0*min(cell(1), cell(2), cell(3))

    if(rcoul.gt.0.0d0)then
      rcoul = min(rcoul, drmin)
    else
      rcoul = drmin
    end if 

    if(rvdw.gt.0.0d0)then
      rvdw = min(rvdw, drmin)
    else
      rvdw = drmin
    end if 

    do i = 1, natom
      do j = 1, 3 
        fa(i, j) = 0.d0
      end do 
      ea(i) = 0.d0 
    end do 

    energy = 0.d0 
    virial = 0.d0 

!      call bonds
    call vdw
!    call coulomb

  end subroutine setForces
!
!- intramolecular bonds
!
!  subroutine bonds
!
!    implicit none
!
!    integer :: i, j, k, l, nx, nk, nl
!    real(8) epot, fr, dx, dy, dz, dr 
!
!
!    nx = 0 
!    do i = 1, nmolecules
!      do j = 1, molecules(i, 1)
!        do k = 1, molecules(i, 2)
!          nk = nx+k
!          do l = k+1, molecules(i, 2)
!            nl = nx+l 
!            dx = ra(nl, 1)-ra(nk, 1)
!            dy = ra(nl, 2)-ra(nk, 2)
!            dz = ra(nl, 3)-ra(nk, 3)
!            call mic(dx, dy, dz)
!            dr = sqrt(dx**2+dy**2+dz**2)
!!            epot = 0.5d0*params(1)*(dr-params(2))**2
!!            fr = -1.0d0*params(1)*(dr-params(2))/dr
!!            fx(nk) = -fr*dx
!!            fy(nk) = -fr*dy
!!            fz(nk) = -fr*dz
!!            ea(nk) = epot
!!            fx(nl) = +fr*dx
!!            fy(nl) = +fr*dy
!!            fz(nl) = +fr*dz
!!            ea(nl) = epot
!          end do 
!        end do 
!        nx = nx+molecules(i, 1)
!      end do 
!    end do
!
!  end subroutine bonds
!
!- Bond constraint rules
!
  subroutine mic(dx, dy, dz)

    implicit none
    real(8), intent(inout) :: dx, dy, dz

    dx = dx - cell(1)*int(2.0d0*dx/cell(1))
    dy = dy - cell(2)*int(2.0d0*dy/cell(2))
    dz = dz - cell(3)*int(2.0d0*dz/cell(3))

  end subroutine mic
!
!- Coulomb interaction
!
!  subroutine coulomb()
!
!    implicit none
!
!    integer :: i, j, nj
!    real(8) :: epot, fr, dx, dy, dz, dr, dvir, rcutoff
!    
!    if(abs(rvdw-rcoul).gt.1.0d-2) then
!      rcutoff = rcoul
!      call neighbour_list(rcutoff)
!    end if
!
!  end subroutine coulomb
!
!- Van der Waals interaction
!
  subroutine vdw()

    implicit none
    integer :: i, k, ni, nj 
    real(8) :: depot, epot, fr, dx, dy, dz, dr
    real(8) :: prm1, prm2, rcutoff
    
    rcutoff = rvdw

    call neighbour_list(rcutoff)

    do i = 1, nlist 
      ni = sites(ilist(i), 1)+1 
      nj = sites(ilist(i), 2)+1
      prm1 = 0.d0
      prm2 = 0.d0
      do k = 1, nvdw 
        if(int(ijinter(k, 1)).eq.atype(ni).and.int(ijinter(k, 2)).eq.atype(nj))then
          prm1 = prmvdw(k, 1)
          prm2 = prmvdw(k, 2)
        end if 
      end do 
      dx = ra(nj, 1)-ra(ni, 1)
      dy = ra(nj, 2)-ra(ni, 2)
      dz = ra(nj, 3)-ra(ni, 3)
      call mic(dx, dy, dz)
      dr = sqrt(dx**2.0d0+dy**2.0d0+dz**2.0d0)
      depot = (prm2/dr)**6.0d0 
      epot = 4.0d0*prm1*(depot-1.0d0)*depot
      fr = 24.d0*prm1*(2.d0*depot-1.d0)*depot/dr**2.0d0 
      fa(ni, 1) = fa(ni, 1)-fr*dx
      fa(ni, 2) = fa(ni, 2)-fr*dy
      fa(ni, 3) = fa(ni, 3)-fr*dz
      fa(nj, 1) = fa(nj, 1)+fr*dx
      fa(nj, 2) = fa(nj, 2)+fr*dy
      fa(nj, 3) = fa(nj, 3)+fr*dz
      ea(ni) = 0.0d0
      virial = virial + fr*dr**2.0d0
      energy = energy + epot
    end do

    call vdw_corr

  end subroutine vdw
!
  subroutine vdw_corr

    implicit none
    integer :: i, j, k, ni, nj
    real(8) :: prm1, prm2, encorr, vircorr, pi 

    pi = acos(-1.d0)

    encorr = 0.d0 
    vircorr = 0.d0 

    j = 1 
    do i = 1, size(nsites)
      prm1 = 0.d0
      prm2 = 0.d0
      ni = sites(j, 1)+1 
      nj = sites(j, 2)+1 
      do k = 1, nvdw
        if(int(ijinter(k, 1)).eq.atype(ni).and.int(ijinter(k, 2)).eq.atype(nj))then
          prm1 = prmvdw(k, 1)
          prm2 = prmvdw(k, 2)
        end if 
      end do 
      encorr = encorr+nsites(i)*4.0d0*prm1*(prm2**12.0d0/(9.0d0*rvdw**9.0d0)-prm2**6.0d0/(3.d0*rvdw**3.0d0))
      vircorr=vircorr+nsites(i)*24.d0*prm1*(2.d0*prm2**12-3.d0*(rvdw*prm2)**6)/(9.d0*rvdw**9)
      j = j+nsites(i)
    end do 
!
    encorr = 2.0d0*pi*encorr/(cell(1)*cell(2)*cell(3))
    vircorr = 2.0d0*pi*vircorr/(cell(1)*cell(2)*cell(3))

    energy = energy + encorr
    virial = virial + vircorr

  end subroutine vdw_corr
!
!-Neighbours list subroutine function
!
  subroutine neighbour_list(rcutoff)

    implicit none
    integer :: i, j, k, nx, ni, nj
    real(8) :: dx, dy, dz, dr
    real(8), intent(in) :: rcutoff

    nx = 1 
    k = 1  
    do i = 1, size(nsites)
      do j = 1, nsites(i)
        ni = sites(k, 1)+1 
        nj = sites(k, 2)+1 
        dx = ra(nj, 1)-ra(ni, 1)
        dy = ra(nj, 2)-ra(ni, 2)
        dz = ra(nj, 3)-ra(ni, 3)
        call mic(dx, dy, dz)
        dr = sqrt(dx**2+dy**2+dz**2)
        if (dr.lt.rcutoff) then
          ilist(nx) = k
          nx = nx+1
        end if 
        k = k+1
      end do
    end do
                 
    nlist = nx-1

  end subroutine neighbour_list
!
!  subroutine vdw()
!
!    implicit none
!    integer :: i, j, k, nj 
!    real(8) :: depot, epot, fr, dx, dy, dz, dr
!    real(8) :: prm1, prm2, rcutoff
!    
!    rcutoff = rvdw
!
!    call neighbour_list(rcutoff)
!
!    do i = 1, size(nlist)
!      do j = 1, nlist(i)
!        nj = ilist(i, j)
!        prm1 = 0.d0
!        prm2 = 0.d0
!        do k = 1, nvdw 
!          if(int(ijinter(k, 1)).eq.atype(i).and.int(ijinter(k, 2)).eq.atype(nj))then
!            prm1 = prmvdw(k, 1)
!            prm2 = prmvdw(k, 2)
!          end if 
!        end do 
!        dx = ra(nj, 1)-ra(i, 1)
!        dy = ra(nj, 2)-ra(i, 2)
!        dz = ra(nj, 3)-ra(i, 3)
!        call mic(dx, dy, dz)
!        dr = sqrt(dx**2.0d0+dy**2.0d0+dz**2.0d0)
!        depot = (prm2/dr)**6.0d0 
!        epot = 4.0d0*prm1*(depot-1.0d0)*depot
!        fr = 24.d0*prm1*(2.d0*depot-1.d0)*depot/dr**2.0d0 
!        fa(i, 1) = fa(i, 1)-fr*dx
!        fa(i, 2) = fa(i, 2)-fr*dy
!        fa(i, 3) = fa(i, 3)-fr*dz
!        fa(nj, 1) = fa(nj, 1)+fr*dx
!        fa(nj, 2) = fa(nj, 2)+fr*dy
!        fa(nj, 3) = fa(nj, 3)+fr*dz
!        ea(i) = ea(i)+epot
!        virial = virial + fr*dr**2.0d0
!      end do
!      energy = energy + ea(i)
!    end do
!
!    call vdw_corr
!
!  end subroutine vdw
!!
!  subroutine vdw_corr
!
!    implicit none
!    integer :: i, j, k, ni, nj
!    real(8) :: prm1, prm2, encorr, vircorr, pi 
!
!    pi = acos(-1.d0)
!
!    encorr = 0.d0 
!    vircorr = 0.d0 
!
!    do i = 1, nsites
!      do j = i, nsites
!        prm1 = 0.d0
!        prm2 = 0.d0
!        do k = 1, nvdw
!          if(int(ijinter(k, 1)).eq.sites(i, 1).and.int(ijinter(k, 2)).eq.sites(j, 1))then
!            prm1 = prmvdw(k, 1)
!            prm2 = prmvdw(k, 2)
!          end if 
!        end do 
!        ni = sites(i, 2)
!        nj = sites(j, 2)
!        encorr = encorr+ni*nj*4.0d0*prm1*(prm2**12.0d0/(9.0d0*rvdw**9.0d0)-prm2**6.0d0/(3.d0*rvdw**3.0d0))
!        vircorr=vircorr+ni*nj*24.d0*prm1*(2.d0*prm2**12-3.d0*(rvdw*prm2)**6)/(9.d0*rvdw**9)
!      end do 
!    end do 
!!
!    encorr = 2.0d0*pi*encorr/(cell(1)*cell(2)*cell(3))
!    vircorr = 2.0d0*pi*vircorr/(cell(1)*cell(2)*cell(3))
!
!    energy = energy + encorr
!    virial = virial + vircorr
!
!  end subroutine vdw_corr
!  subroutine neighbour_list(rcutoff)
!
!    implicit none
!    integer :: i, j, nx, i1, i2, j1, j2, g, gg, p, pp, ix 
!    real(8) :: dx, dy, dz, dr
!    real(8), intent(in) :: rcutoff
!
!    if (allocated(nlist).eqv..FALSE.) then
!      allocate(nlist(natom), ilist(natom, natom))
!    end if 
!
!    gg = 0
!    ix = 0 
!    do i = 1, nmolecules
!      do i1 = 1, molecules(i, 1)
!        do i2 = 1, molecules(i, 2)
!          nx = 1
!          ix = ix+1 
!          g = gg+i2+(i1-1)*molecules(i, 2)
!          do j1 = i1+1, molecules(i, 1)
!            do j2 = 1, molecules(i, 2)
!              p = gg+j2+(j1-1)*molecules(i, 2)
!              dx = ra(p, 1)-ra(g, 1)
!              dy = ra(p, 2)-ra(g, 2)
!              dz = ra(p, 3)-ra(g, 3)
!              call mic(dx, dy, dz)
!              dr = sqrt(dx**2+dy**2+dz**2)
!              if (dr.lt.rcutoff) then
!                ilist(g, nx) = p
!                nx = nx+1
!              end if 
!            end do
!          end do 
!          nlist(g) = nx-1
!        end do 
!      end do 
!      gg = gg+molecules(i, 1)*molecules(i, 2)
!    end do 
!    gg = 0 
!    ix = 0
!    do i = 1, nmolecules
!      do i1 = 1, molecules(i, 1)
!        do i2 = 1, molecules(i, 2)
!          nx = 1
!          ix = ix+1 
!          g = gg+i2+(i1-1)*molecules(i, 2)
!          pp = gg+molecules(i, 1)*molecules(i, 2)
!          do j = i+1, nmolecules
!            do j1 = 1, molecules(j, 1)
!              do j2 = 1, molecules(j, 2)
!                p = pp+j2+(j1-1)*molecules(j, 2)
!                dx = ra(p, 1)-ra(g, 1)
!                dy = ra(p, 2)-ra(g, 2)
!                dz = ra(p, 3)-ra(g, 3)
!                call mic(dx, dy, dz)
!                dr = sqrt(dx**2+dy**2+dz**2)
!                if (dr.lt.rcutoff)then
!                  ilist(g, nx) = p
!                  nx = nx+1
!                end if 
!              end do
!            end do
!            pp = pp+molecules(j, 1)*molecules(j, 2)
!          end do
!          nlist(g) = nlist(g)+nx-1
!        end do
!      end do
!      gg = gg+molecules(i, 1)*molecules(i, 2)
!    end do
!                 
!  end subroutine neighbour_list

end module force_field_module
