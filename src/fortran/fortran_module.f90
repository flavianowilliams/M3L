module fortran_module 
  integer :: natom, nmolecules, nsites, nvdw, nfree
  integer, allocatable, dimension(:) :: nlist
  integer, allocatable, dimension(:,:) :: ilist, molecules, sites 
  real(8) :: timestep, tstat, pstat, bfactor, temp_bath, press_bath, friction
  real(8) :: ekinetic, energy, temperature, pressure, virial, rvdw, rcoul, temp_friction, press_friction, volume 
  real(8), dimension(3) :: cell
  real(8), allocatable, dimension(:, :) :: params, atom
contains
!
!- preparing variables
!
  subroutine prepare

    real(8) :: drmin

    drmin = 0.5d0*min(cell(1), cell(2), cell(3))

    if(rvdw.gt.0.d0)then
      rvdw = min(rvdw, drmin)
    else
      rvdw = drmin
    end if 

    if (allocated(params).eqv..FALSE.) then
      allocate(params(nvdw+1, 4))
    end if 

  end subroutine prepare
!
!-Calculate forces a atomic energy
!
  subroutine forces

    implicit none
    integer :: i

    do i = 1, natom
      atom(i, 10) = 0.d0
      atom(i, 11) = 0.d0
      atom(i, 12) = 0.d0
      atom(i, 13) = 0.d0
    end do 

!      call bonds
    call vdw

  end subroutine forces
!
!- intramolecular bonds
!
!  subroutine bonds
!    implicit none
!    real(8) epot, fr, dx, dy, dz, dr 
!    
!    dx = rx(2)-rx(1)
!    dy = ry(2)-ry(1)
!    dz = rz(2)-rz(1)
!    call mic(dx, dy, dz)
!    dr = sqrt(dx**2+dy**2+dz**2)
!    epot = 0.5d0*params(1)*(dr-params(2))**2
!    fr = -1.0d0*params(1)*(dr-params(2))/dr
!    fx(1) = -fr*dx
!    fy(1) = -fr*dy
!    fz(1) = -fr*dz
!    ea(1) = epot
!    fx(2) = +fr*dx
!    fy(2) = +fr*dy
!    fz(2) = +fr*dz
!    ea(2) = epot
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
!- Van der Waals interaction
!
  subroutine vdw

    implicit none
    integer :: i, j, k, nj
    real(8) :: depot, epot, fr, dx, dy, dz, dr, dvir
    real(8) :: prm1, prm2, p1a, p1b, p2a, p2b
    
    call neighbour_list

    energy = 0.d0 
    virial = 0.d0 
    do i = 1, size(nlist)
      dvir = 0.d0 
      atom(i, 13) = 0.d0
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
        dvir = dvir+fr*dr**2.0d0 
      end do
      energy = energy + atom(i, 13)
      virial = virial + dvir
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
  subroutine neighbour_list

    implicit none
    integer :: i, j, nx, i1, i2, j1, j2, g, gg, p, pp, ix 
    real(8) :: dx, dy, dz, dr

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
              if (dr.lt.rvdw) then
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
                if (dr.lt.rvdw) then
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
!
!-Update volume 
!
  subroutine setVolume()

    volume = cell(1)*cell(2)*cell(3)    

  end subroutine
!
!-Ensemble nvt berendsen
!
  subroutine ekinetic_func()

    implicit none
    integer :: i 
    real(8) :: sum

    sum = 0.d0
    do i = 1, natom
      sum = sum+atom(i, 2)*(atom(i, 7)**2+atom(i, 8)**2+atom(i, 9)**2)
    end do

    ekinetic = 0.5d0*sum

  end subroutine ekinetic_func

  subroutine temperature_func
  
    temperature = 2.0d0*ekinetic/nfree

  end subroutine temperature_func

  subroutine pressure_func

    implicit none 
    
    pressure = (2.0d0*ekinetic+virial)/(3.0d0*volume)

  end subroutine pressure_func
!
!-Ensemble nvt berendsen
!
!  subroutine nvt_berendsen
!
!    implicit none
!    integer :: i 
!    real(8) :: qui, sigma
!    
!    do i = 1, natom 
!      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
!      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
!      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
!      rx(i) = rx(i)+vx(i)*timestep
!      ry(i) = ry(i)+vy(i)*timestep
!      rz(i) = rz(i)+vz(i)*timestep
!
!    end do 
!
!    call ccp
!
!    call forces
!
!    do i = 1, natom
!      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
!      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
!      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
!    end do 
!
!    call ekinetic_func
!
!    sigma = 0.5d0*nfree*temp_bath
!    qui = sqrt(1.d0+timestep*(sigma/ekinetic-1.0d0)/tstat)
!
!    do i = 1, natom
!      vx(i) = vx(i)*qui
!      vy(i) = vy(i)*qui
!      vz(i) = vz(i)*qui
!    end do 
!
!    call ekinetic_func
!    call temperature_func
!
!  end subroutine nvt_berendsen

  subroutine friction_func

    implicit none

    real(8) :: sigma

    sigma = 0.5d0*nfree*temp_bath
    friction = friction+0.25d0*timestep*(ekinetic-sigma)/(sigma*tstat**2)

  end subroutine friction_func

  subroutine setTempFriction

    implicit none

    real(8) :: term_mass, bar_mass, func

    term_mass = nfree*temp_bath*tstat**2.d0 
    bar_mass = nfree*temp_bath*pstat**2.0d0

    func = 0.5d0*timestep*(nfree*(temperature-temp_bath)+(bar_mass*press_friction**2.0d0-temp_bath))/term_mass

    temp_friction = temp_friction+func

  end subroutine setTempFriction

  subroutine press_friction_func

    implicit none

    real(8) :: term_mass, bar_mass, func

    term_mass = nfree*temp_bath*tstat**2.d0 
    bar_mass = nfree*temp_bath*pstat**2.0d0

    func = 0.5d0*timestep*(3.0d0*volume*(pressure-press_bath)/bar_mass-temp_friction*press_friction)

    press_friction = press_friction + func

  end subroutine press_friction_func

  subroutine eta_func(eta)

    implicit none
    real(8), intent(out) :: eta

    eta = (1.0d0-bfactor*timestep*(press_bath-pressure)/pstat)**(1.0d0/3.0d0)

  end subroutine eta_func

!  subroutine nvt_hoover
!
!    implicit none
!    integer :: i 
!    
!    call ekinetic_func
!    call friction_func
!
!    do i = 1, natom 
!      vx(i) = vx(i)*exp(-0.5d0*timestep*friction)
!      vy(i) = vy(i)*exp(-0.5d0*timestep*friction)
!      vz(i) = vz(i)*exp(-0.5d0*timestep*friction)
!    end do 
!
!    call ekinetic_func
!    call friction_func
!
!    do i = 1, natom
!      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
!      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
!      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
!      rx(i) = rx(i)+vx(i)*timestep
!      ry(i) = ry(i)+vy(i)*timestep
!      rz(i) = rz(i)+vz(i)*timestep
!    end do 
!
!    call ccp
!    call forces
!
!    do i = 1, natom
!      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
!      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
!      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
!    end do 
!
!    call ekinetic_func
!    call friction_func
!
!    do i = 1, natom
!      vx(i) = vx(i)*exp(-0.5d0*timestep*friction)
!      vy(i) = vy(i)*exp(-0.5d0*timestep*friction)
!      vz(i) = vz(i)*exp(-0.5d0*timestep*friction)
!    end do 
!
!    call ekinetic_func
!    call friction_func
!    call temperature_func
!    call pressure_func
!
!  end subroutine nvt_hoover

  subroutine npt_berendsen

    implicit none
    integer :: i
    real(8) :: eta, qui, sigma

    call eta_func(eta)

    do i = 1, natom 
      atom(i, 7) = atom(i, 7)+atom(i, 10)*0.5d0*timestep/atom(i, 2)
      atom(i, 8) = atom(i, 8)+atom(i, 11)*0.5d0*timestep/atom(i, 2)
      atom(i, 9) = atom(i, 9)+atom(i, 12)*0.5d0*timestep/atom(i, 2)
      atom(i, 4) = atom(i, 4)*eta+atom(i, 7)*timestep
      atom(i, 5) = atom(i, 5)*eta+atom(i, 8)*timestep
      atom(i, 6) = atom(i, 6)*eta+atom(i, 9)*timestep
    end do 

    do i = 1,3
      cell(i) = cell(i)*eta 
    end do 

    call setVolume
    call ccp
    call forces

    do i = 1, natom 
      atom(i, 7) = atom(i, 7)+atom(i, 10)*0.5d0*timestep/atom(i, 2)
      atom(i, 8) = atom(i, 8)+atom(i, 11)*0.5d0*timestep/atom(i, 2)
      atom(i, 9) = atom(i, 9)+atom(i, 12)*0.5d0*timestep/atom(i, 2)
    end do 

    call ekinetic_func

    sigma = 0.5d0*nfree*temp_bath
    qui = sqrt(1.0d0+timestep*(sigma/ekinetic-1.0d0)/tstat)

    do i = 1, natom 
      atom(i, 7) = atom(i, 7)*qui
      atom(i, 8) = atom(i, 8)*qui
      atom(i, 9) = atom(i, 9)*qui
    end do 

    call ekinetic_func
    call temperature_func
    call pressure_func

  end subroutine npt_berendsen

!  subroutine npt_hoover
!
!    implicit none
!
!    integer :: i
!
!    call setTempFriction
!
!    do i = 1, natom 
!      vx(i) = vx(i)-0.5d0*timestep*temp_friction*vx(i)
!      vy(i) = vy(i)-0.5d0*timestep*temp_friction*vy(i)
!      vz(i) = vz(i)-0.5d0*timestep*temp_friction*vz(i)
!    end do 
!
!    call press_friction_func
!    
!    do i = 1, natom 
!      vx(i) = vx(i)-0.5d0*timestep*temp_friction*vx(i)
!      vy(i) = vy(i)-0.5d0*timestep*temp_friction*vy(i)
!      vz(i) = vz(i)-0.5d0*timestep*temp_friction*vz(i)
!    end do 
!
!    do i = 1, natom 
!      vx(i) = vx(i)+fx(i)*0.5d0*timestep/mass(i)
!      vy(i) = vy(i)+fy(i)*0.5d0*timestep/mass(i)
!      vz(i) = vz(i)+fz(i)*0.5d0*timestep/mass(i)
!    end do 
!
!    do i = 1, natom 
!      rx(i) = rx(i)+vx(i)*timestep
!      ry(i) = ry(i)+vy(i)*timestep
!      rz(i) = rz(i)+vz(i)*timestep
!    end do 
!
!  end subroutine npt_hoover

end module fortran_module
