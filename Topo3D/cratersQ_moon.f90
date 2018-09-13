program cratersQ_Moon
!****************************************************************************
!   cratersQ: program to calculate insolation with shadowing and reflections
!             for airless body
!  
!   written by Norbert Schorghofer 2010-2016 
!****************************************************************************
  use filemanager
  use allinterfaces
  use newhorizons
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8  
  real(8), parameter :: solarDay = 29.53*86400.

  integer nsteps, n, i, j, nm, k, CCMAX, iii, jjj, nm2
  real(8) tmax, dt, latitude, dtsec
  real(8) R, Decl, HA, sdays
  real(8) azSun, sinbeta, smax, emiss, v
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qn, QIR, Qrefl   ! incoming
  integer, dimension(NSx,NSy) :: cc
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: dO12  
  real(8), dimension(NSx,NSy) :: Tsurf, Qvis, skysize, Qabs, albedo, QIRin, QIRre
  real(8) Qmeans(NSx,NSy,4)
  real(8), dimension(NSx,NSy) :: Qmax1, Qmax2, Tmean, Tmaxi, Tb, Tmaxi2
  real(8), allocatable :: T(:,:,:), Qnm1(:,:)  ! subsurface
  logical, parameter :: reflection=.true., subsurface=.false.
  
  dt=0.01d0 
  tmax = 2.
  latitude = 87.
  ! azimuth in degrees east of north, 0=north facing
  albedo(:,:) = 0.12d0
  emiss = 0.95d0

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'Reflections:',reflection,'Subsurface:',subsurface

  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r
  Tsurf=0.; Qrefl=0.; QIRre=0.  
  Qmeans=0.; Tmean=0.; nm=0
  Qmax1=0.; Qmax2=0.; 
  Tmaxi=0.; Tmaxi2=0.
  Tb=0.; nm2=0
  
  print *,'...reading horizons file...'
  call readhorizons('horizons.'//fileext)

  if (reflection) then
     print *,'...reading huge fieldofviews file...'
     call getmaxfieldsize(NSx,NSy,fileext,CCMAX)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), dO12(NSx,NSy,CCMAX))
     call getfieldofview(NSx,NSy,fileext,cc,ii,jj,dO12,skysize,CCMAX)
  endif
  if (subsurface) allocate(T(1000,NSx,NSy), Qnm1(NSx,NSy))

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     !print *,sdays
     Decl=0.; R=1.

     ! for small landscapes these can be here
     HA=2.*pi*mod(sdays,1.)   ! hour angle
     call equatorial2horizontal(Decl,latitude,HA,sinbeta,azSun)

     print *,sdays,n,HA/d2r,asin(sinbeta)/d2r,azSun/d2r

     do i=2,NSx-1
        do j=2,NSy-1
           smax = getonehorizon(i,j,azSun)
           Qn(i,j)=flux_wshad(R,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),smax)
        enddo
     enddo

     if (reflection) then
        Qvis(:,:) = Qn + Qrefl
        QIRin(:,:) = QIR + QIRre
        do i=2,NSx-1
           do j=2,NSy-1
              QIR(i,j)=0.; Qrefl(i,j)=0.; QIRre(i,j)=0.
              do k=1,cc(i,j)
                 iii = ii(i,j,k); jjj = jj(i,j,k)
                 v = viewing_angle(i,j,iii,jjj,h)
                 Qrefl(i,j) = Qrefl(i,j) + dO12(i,j,k)/pi*cos(v)*albedo(iii,jjj)*Qvis(iii,jjj)
                 QIR(i,j) = QIR(i,j) + dO12(i,j,k)/pi*cos(v)*emiss*sigSB*Tsurf(iii,jjj)**4
                 QIRre(i,j) = QIRre(i,j) + dO12(i,j,k)/pi*cos(v)*(1-emiss)*QIRin(iii,jjj)
              enddo
           enddo
        enddo
        Qabs(:,:)=(1.-albedo(:,:))*(Qn+Qrefl)+emiss*(QIR+QIRre)  ! Q absorbed
     else
        Qabs(:,:)=(1.-albedo(:,:))*Qn
     endif

     if (subsurface) then
        if (n==0) then ! initialization
           Qnm1 = Qabs
           !Tsurf = 0.
           Tsurf = -200. ! negative of initialization temperature
        endif

        ! for faster convergence
        if (n>=10 .and. n<10+nint(1/dt)) then
           Tb = Tb+Tsurf
           nm2 = nm2+1
        endif
        if (n==10+nint(1/dt)) then
           Tb = Tb/nm2
           Tsurf = -Tb  ! re-initialize
        endif

        ! main part
        do i=2,NSx-1
           do j=2,NSy-1
              call subsurfaceconduction(T(:,i,j),Tsurf(i,j),dtsec,Qnm1(i,j),Qabs(i,j),emiss,solarDay)
           enddo
        enddo
        Qnm1 = Qabs
     else
        Tsurf = (Qabs/sigSB/emiss)**0.25
     endif
     
     if (sdays > tmax-1.) then
        Qmeans(:,:,1) = Qmeans(:,:,1) + Qn
        where (Qn>Qmax1) Qmax1=Qn  ! maximum direct
        Qmeans(:,:,2) = Qmeans(:,:,2) + Qabs
        where (Qabs>Qmax2) Qmax2=Qabs  ! maximum total
        Qmeans(:,:,3) = Qmeans(:,:,3) + QIR
        Qmeans(:,:,4) = Qmeans(:,:,4) + Qrefl
        Tmean = Tmean + Tsurf
        where (Tsurf>Tmaxi)
           Tmaxi2=Tmaxi
           Tmaxi=Tsurf
        elsewhere (Tsurf>Tmaxi2)
           Tmaxi2=Tsurf
        end where
        nm=nm+1
     endif

  enddo  ! end of time loop

  if (reflection) deallocate(ii, jj, dO12)
  if (subsurface) deallocate(T, Qnm1)
  
  Qmeans=Qmeans/nm
  Tmean=Tmean/nm

  open(unit=21,file='qmean.dat',status='unknown',action='write')
  open(unit=22,file='qinst.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,4(1x,f6.1),1x,f5.1,2(1x,f6.1),1x,f5.1)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmeans(i,j,1:4), &
             & Tmean(i,j),Qmax1(i,j),Qmax2(i,j),Tmaxi(i,j)
        write(22,'(2(i4,1x),f9.2,2x,f6.4,4(1x,f6.1),1x,f5.1)') &  ! instanteneous values
             & i,j,h(i,j),surfaceSlope(i,j),Qn(i,j),Qabs(i,j),Qir(i,j),Qrefl(i,j), &
             & Tsurf(i,j)
     enddo
  enddo
  close(21)
  close(22)
end program cratersQ_Moon
 


subroutine subsurfaceconduction(T,Tsurf,dtsec,Qn,Qnp1,emiss,solarDay)
  ! 1d subsurface conduction
  use allinterfaces, only : conductionQ
  implicit none
  integer, parameter :: NMAX=1000, Ni=5
  real(8), parameter :: pi=3.1415926535897932
  real(8), parameter :: Fgeotherm = 0.029
  integer, parameter :: nz=25
  real(8), intent(INOUT) :: T(NMAX), Tsurf
  real(8), intent(IN) :: dtsec,Qn,Qnp1,emiss,solarDay
  integer i, k
  real(8) zmax, zfac, Fsurf, Tinit, delta
  real(8) Tsurfold, Told(1:nz), Qarti, Qartiold 
  real(8), save :: ti(NMAX), rhocv(NMAX), z(NMAX)
  logical, save :: first = .true.

  if (first) then ! initialize grid
     ti(:) = 100.;  rhocv(:) = 1200.*800.  ! adjust
     zmax=0.5; zfac = 1.05  ! adjust

     delta = ti(1)/rhocv(1)*sqrt(solarDay/pi)  ! skin depth

     call setgrid(nz,z,zmax,zfac)
     if (z(6)>delta) then
        print *,'WARNING: less than 6 points within diurnal skin depth'
     endif
     do i=1,nz
        if (z(i)<delta) cycle
        print *,i-1,' grid points within diurnal skin depth'
        exit
     enddo
     if (z(1)<1.e-5) print *,'WARNING: first grid point is too shallow'
     open(unit=30,file='z',status='unknown');
     write(30,*) (z(i),i=1,nz)
     close(30)

     write(*,*) 'Subsurface model parameters'
     write(*,*) '   nz=',nz,' zmax=',zmax,' zfac=',zfac
     write(*,*) '   Thermal inertia=',ti(1),' rho*c=',rhocv(1)
     print *,'   Diurnal skin depth=',delta,' Geothermal flux=',Fgeotherm

     first = .false.
  endif
  
  if (Tsurf<=0.) then  ! initialize temperature profile
     if (Tsurf==0.) Tinit=273.
     if (Tsurf<0.) Tinit=-Tsurf
     T(1:nz) = Tinit
     Tsurf = Tinit
  endif

  Tsurfold=Tsurf
  Told(1:nz)=T(1:nz)
  call conductionQ(nz,z,dtsec,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeotherm,Fsurf)
  
  ! fixes rapid sunrise on sloped surface
  ! artificial flux smoothing
  if (Tsurf>2*Tsurfold .or. Tsurf<Tsurfold/2) then
     Tsurf = Tsurfold; T(1:nz) = Told(1:nz)
     do k=1,Ni
        Qartiold = ((Ni-k+1)*Qn + (k-1)*Qnp1)/real(Ni)
        Qarti = ((Ni-k)*Qn + k*Qnp1)/real(Ni)
        call conductionQ(nz,z,dtsec/Ni,Qartiold,Qarti,T,ti,rhocv,emiss,Tsurf,Fgeotherm,Fsurf)
     enddo
  endif
  
end subroutine subsurfaceconduction
