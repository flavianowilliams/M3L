module force_field_module 
  integer :: natom, nmolecules, nsites, nvdw
  integer, allocatable, dimension(:) :: nlist
  integer, allocatable, dimension(:,:) :: ilist, molecules, sites 
  real(8) :: virial, rvdw, rcoul, energy 
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: params, atom
contains
!
!-Calculate forces a atomic energy
!
  subroutine setForces()

    implicit none
    integer :: i
    real(8) :: drmin

    drmin = 0.5d0*min(cell(1), cell(2), cell(3))

    if(rcoul.gt.0.d0)then
      rcoul = min(rcoul, drmin)
    else
      rcoul = drmin
    end if 

    if(rvdw.gt.0.d0)then
      rvdw = min(rvdw, drmin)
    else
      rvdw = drmin
    end if 

    if (allocated(params).eqv..FALSE.) then
      allocate(params(nvdw+1, 4))
    end if 

    do i = 1, natom
      atom(i, 10) = 0.d0
      atom(i, 11) = 0.d0
      atom(i, 12) = 0.d0
      atom(i, 13) = 0.d0
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
  subroutine bonds

    implicit none

    integer :: i, j, k, l, nx, nk, nl
    real(8) epot, fr, dx, dy, dz, dr 


    nx = 0 
    do i = 1, nmolecules
      do j = 1, molecules(i, 1)
        do k = 1, molecules(i, 2)
          nk = nx+k
          do l = k+1, molecules(i, 2)
            nl = nx+l 
            dx = atom(nl, 4)-atom(nk, 4)
            dy = atom(nl, 5)-atom(nk, 5)
            dz = atom(nl, 6)-atom(nk, 6)
            call mic(dx, dy, dz)
            dr = sqrt(dx**2+dy**2+dz**2)
!            epot = 0.5d0*params(1)*(dr-params(2))**2
!            fr = -1.0d0*params(1)*(dr-params(2))/dr
!            fx(nk) = -fr*dx
!            fy(nk) = -fr*dy
!            fz(nk) = -fr*dz
!            ea(nk) = epot
!            fx(nl) = +fr*dx
!            fy(nl) = +fr*dy
!            fz(nl) = +fr*dz
!            ea(nl) = epot
          end do 
        end do 
        nx = nx+molecules(i, 1)
      end do 
    end do

  end subroutine bonds
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
  subroutine coulomb()

    implicit none

    integer :: i, j, nj
    real(8) :: epot, fr, dx, dy, dz, dr, dvir, rcutoff
    
    if(abs(rvdw-rcoul).gt.1.0d-2) then
      rcutoff = rcoul
      call neighbour_list(rcutoff)
    end if

    do i = 1, size(nlist)
      dvir = 0.d0 
      atom(i, 13) = 0.d0
      do j = 1, nlist(i)
        nj = ilist(i, j)
        dx = atom(nj, 4)-atom(i, 4)
        dy = atom(nj, 5)-atom(i, 5)
        dz = atom(nj, 6)-atom(i, 6)
        call mic(dx, dy, dz)
        dr = sqrt(dx**2.0d0+dy**2.0d0+dz**2.0d0)
        epot = atom(i, 3)*atom(nj, 3)*(1.d0/dr+dr/rcoul**2.0d0-2.0d0/rcoul)
        fr = atom(i, 3)*atom(nj, 3)*(1.0d0/dr**2.0d0-1.0d0/rcoul**2.0d0)/dr
        atom(i, 10) = atom(i, 10)-fr*dx
        atom(i, 11) = atom(i, 11)-fr*dy
        atom(i, 12) = atom(i, 12)-fr*dz
        atom(nj, 10) = atom(nj, 10)+fr*dx
        atom(nj, 11) = atom(nj, 11)+fr*dy
        atom(nj, 12) = atom(nj, 12)+fr*dz
        atom(i, 13) = atom(i, 13)+epot
        dvir = dvir+fr*dr**2.0d0 
      end do
      energy = energy + atom(i, 13)
      virial = virial + dvir
    end do

  end subroutine coulomb
!
!- Van der Waals interaction
!
  subroutine vdw()

    implicit none
    integer :: i, j, k, nj 
    real(8) :: depot, epot, fr, dx, dy, dz, dr
    real(8) :: prm1, prm2, p1a, p1b, p2a, p2b, rcutoff
    
    rcutoff = rvdw

    call neighbour_list(rcutoff)

!    dvir = 0.d0 
    do i = 1, size(nlist)
!      atom(i, 13) = 0.d0
      p1a = 0.d0
      p2a = 0.d0 
      do k = 1, nvdw 
        if(int(params(k+1, 1)).eq.int(atom(i, 1)))then
          p1a = params(k+1, 3)
          p2a = params(k+1, 4)
        end if 
      end do 
      do j = 1, nlist(i)
        nj = ilist(i, j)
        p1b = 0.d0 
        p2b = 0.d0
        do k = 1, nvdw 
          if(int(params(k+1, 1)).eq.int(atom(nj, 1)))then
            p1b = params(k+1, 3)
            p2b = params(k+1, 4)
          end if 
        end do 
        prm1 = sqrt(p1a*p1b)
        prm2 = 0.5d0*(p2a+p2b)
        dx = atom(nj, 4)-atom(i, 4)
        dy = atom(nj, 5)-atom(i, 5)
        dz = atom(nj, 6)-atom(i, 6)
        call mic(dx, dy, dz)
        dr = sqrt(dx**2.0d0+dy**2.0d0+dz**2.0d0)
        depot = (prm2/dr)**6.0d0 
        epot = 4.0d0*prm1*(depot-1.0d0)*depot
        fr = 24.d0*prm1*(2.d0*depot-1.d0)*depot/dr**2.0d0 
        atom(i, 10) = atom(i, 10)-fr*dx
        atom(i, 11) = atom(i, 11)-fr*dy
        atom(i, 12) = atom(i, 12)-fr*dz
        atom(nj, 10) = atom(nj, 10)+fr*dx
        atom(nj, 11) = atom(nj, 11)+fr*dy
        atom(nj, 12) = atom(nj, 12)+fr*dz
        atom(i, 13) = atom(i, 13)+epot
        virial = virial + fr*dr**2.0d0
!        dvir = dvir+fr*dr**2.0d0 
      end do
      energy = energy + atom(i, 13)
!      virial = virial + dvir
    end do

    call vdw_corr

  end subroutine vdw
!
  subroutine vdw_corr

    implicit none
    integer :: i, j, k, ni, nj
    real(8) :: prm1, prm2, prm1a, prm1b, prm2a, prm2b, encorr, vircorr, pi 

    pi = acos(-1.d0)

    encorr = 0.d0 
    vircorr = 0.d0 

    do i = 1, nsites
      prm1a = 0.d0 
      prm2a = 0.d0 
      do k = 1, nvdw 
        if(params(k+1, 1).eq.sites(i, 1))then
          prm1a = params(k+1, 3)
          prm2a = params(k+1, 4)
        end if 
      end do 
      do j = i, nsites
        prm1b = 0.d0 
        prm2b = 0.d0
        do k = 1, nvdw
          if(params(k+1, 1).eq.sites(j, 1))then
            prm1b = params(k+1, 3)
            prm2b = params(k+1, 4)
          end if 
        end do 
        prm1 = sqrt(prm1a*prm1b)
        prm2 = 0.5d0*(prm2a+prm2b)
        ni = sites(i, 2)
        nj = sites(j, 2)
        encorr = encorr+ni*nj*4.0d0*prm1*(prm2**12.0d0/(9.0d0*rvdw**9.0d0)-prm2**6.0d0/(3.d0*rvdw**3.0d0))
        vircorr = vircorr-ni*nj*24.0d0*prm1*(2.0d0*prm2**12.0d0/(9.0d0*rvdw**9.0d0)-prm2**6.0d0/(3.d0*rvdw**3.0d0))
      end do 
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
    integer :: i, j, nx, i1, i2, j1, j2, g, gg, p, pp, ix 
    real(8) :: dx, dy, dz, dr
    real(8), intent(in) :: rcutoff

    if (allocated(nlist).eqv..FALSE.) then
      allocate(nlist(natom), ilist(natom, natom))
    end if 

    gg = 0
    ix = 0 
    do i = 1, nmolecules
      do i1 = 1, molecules(i, 1)
        do i2 = 1, molecules(i, 2)
          nx = 1
          ix = ix+1 
          g = gg+i2+(i1-1)*molecules(i, 2)
          do j1 = i1+1, molecules(i, 1)
            do j2 = 1, molecules(i, 2)
              p = gg+j2+(j1-1)*molecules(i, 2)
              dx = atom(p, 4)-atom(g, 4)
              dy = atom(p, 5)-atom(g, 5)
              dz = atom(p, 6)-atom(g, 6)
              call mic(dx, dy, dz)
              dr = sqrt(dx**2+dy**2+dz**2)
              if (dr.lt.rcutoff) then
                ilist(g, nx) = p
                nx = nx+1
              end if 
            end do
          end do 
          nlist(g) = nx-1
        end do 
      end do 
      gg = gg+molecules(i, 1)*molecules(i, 2)
    end do 
    gg = 0 
    ix = 0
    do i = 1, nmolecules
      do i1 = 1, molecules(i, 1)
        do i2 = 1, molecules(i, 2)
          nx = 1
          ix = ix+1 
          g = gg+i2+(i1-1)*molecules(i, 2)
          pp = gg+molecules(i, 1)*molecules(i, 2)
          do j = i+1, nmolecules
            do j1 = 1, molecules(j, 1)
              do j2 = 1, molecules(j, 2)
                p = pp+j2+(j1-1)*molecules(j, 2)
                dx = atom(p, 4)-atom(g, 4)
                dy = atom(p, 5)-atom(g, 5)
                dz = atom(p, 6)-atom(g, 6)
                call mic(dx, dy, dz)
                dr = sqrt(dx**2+dy**2+dz**2)
                if (dr.lt.rcutoff)then
                  ilist(g, nx) = p
                  nx = nx+1
                end if 
              end do
            end do
            pp = pp+molecules(j, 1)*molecules(j, 2)
          end do
          nlist(g) = nlist(g)+nx-1
        end do
      end do
      gg = gg+molecules(i, 1)*molecules(i, 2)
    end do
                 
  end subroutine neighbour_list

end module force_field_module
