module globalparams
    implicit none
    real(8), parameter :: sigSB = 5.6704d-8
    real(8), parameter :: rockdensity = 2500.
    real(8), parameter :: Fgeotherm = 0.018
    integer, parameter :: STEPSPERSOL = 96
end module globalparams

subroutine thermal1d(latdeg,albedo,H,rho_d,slopedeg,azdeg,horEdeg,horWdeg)
!***********************************************************************
! run one lunar thermal model
!***********************************************************************
    use globalparams
    implicit none
    real(8), intent(IN) :: latdeg, albedo, H, rho_d, slopedeg, azdeg, horEdeg, horWdeg
    real(8) Tsurf_out(STEPSPERSOL)
    real(8), parameter :: pi = 3.1415926535897932, d2r = pi/180.
    real(8), parameter :: lunarDay = 86400.*29.53
    integer nz
    real(8) zmax, zfac
    parameter(nz=30, zmax=0.5, zfac=1.07)
    !parameter(nz=40, zmax=0.5, zfac=1.05)
    integer nsteps, n, i, j
    real(8) lat, slope, az, horE, horW
    real(8) rho_s, k_s, k_d, emiss
    real(8) T(nz), tmax, time, Tsurf, Qn, ltime
    real(8) Qnp1, dtsec, sunR, decl, HA
    real(8), dimension(nz) :: ti, rho, cp, z, k
    real(8) Fsurf
    real(8), external :: flux_moon, heatcapacity, conductivity
    logical :: firsttime = .TRUE.
    
    dtsec = lunarDay/STEPSPERSOL
    tmax = 50.  ! lunations
    
    decl = 0.
    !decl = -1.5
    sunR = 1.
    
    print *,'Latitude=',lat,'albedo=',albedo
    if (lat<-90. .or. lat>+90.) stop 'latitude out of range'
    
    nsteps = nint(tmax*lunarDay/dtsec) + STEPSPERSOL/2  ! calculate total number of timesteps

    call setgrid(nz,z,zmax,zfac)
    if (firsttime) then
        open(unit=20,file='z',action='write')
        write(20,'(999(f8.6,1x))') z(1:nz)
        close(20)

        open(unit=21,file='ltimes',action='write')
    end if
    
    ! constant values
    emiss = 0.95
    rho_s = 1100.
    k_s = 7.4e-4; k_d = 3.4e-3 ! k_s = 8.0e-4; k_d = 3.8e-3 ! Martinez et al. (2021) lunar1Dheat values

    block  ! initialize temperatures
        real(8) geof, Tinit
        real(8) Tnom, delta, knz, thIn, ztell

        geof = cos(lat*d2r)/pi
        Tinit = (1370*(1.-albedo)*geof/sigSB)**0.25 - 25. 
        if (geof<=0.) Tinit=50.  ! the poles
        T(1:nz) = Tinit
        Tsurf = Tinit

        ! extra checks
        Tnom = Tinit
        ztell = z(10)
        knz = conductivity(rho_s, Tnom)
        thIn = sqrt( knz*rho_s*heatcapacity(Tnom) )
        delta = (thIn/(rho_s*heatcapacity(Tnom)))*sqrt(lunarDay/pi)  ! skin depth
        
        do i=1,nz
            if (z(i)<delta) cycle
            !print *,i-1,' grid points within diurnal skin depth',delta
            exit
        enddo
    end block

    ! initialize bulk density profile
    do i=1,nz
        rho(i) = rho_d - (rho_d - rho_s) * exp(-z(i)/H)  ! 1100 ... 1800
    enddo

    
    ! time=0, n=0
    HA = 0.
    Qn = flux_moon(sunR,decl*d2r,lat*d2r,HA,albedo)
    j = 0

    do n=0,nsteps-1  ! loop over time steps
        time = (n+1)*dtsec   ! time at n+1
        
        HA = time/lunarDay*2*pi  ! hour angle
        Qnp1 = flux_moon(sunR,decl*d2r,lat*d2r,HA,albedo)

        do i=1,nz ! update thermal properties
            k(i) = conductivity(rho(i), T(i))
            cp(i) = heatcapacity(T(i))
            ti(i) = sqrt( k(i)*rho(i)*cp(i) )
        enddo
        
        call conductionQ(nz,z,dtsec,Qn,Qnp1,T(:),ti(:),rho(:)*cp(:),emiss,Tsurf,Fgeotherm,Fsurf)
        Qn = Qnp1
        
        if (n >= nsteps-STEPSPERSOL) then  ! last lunation
        ltime = mod(time/lunarDay*360-1e-10 + 180, 360.d0)
        !write(22,'(f6.3,1x,f8.3,1(1x,f7.2))') ltime*24./360.,Qn,Tsurf
        j = j+1
        Tsurf_out(j) = Tsurf
        if (firsttime) write(21,'(f6.3)') ltime*24./360.
        end if
    end do

    if (firsttime) close(21)
    write(22,'(*(1x,f7.2))') Tsurf_out

    firsttime = .FALSE.
end subroutine thermal1d
    
    