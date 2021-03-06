!***********************************************************************
! a collection of commonly-used functions
!***********************************************************************

function flux2T(Q,albedo,emiss)
  implicit none
  real(8), intent(IN) :: Q, emiss, albedo
  real*8, parameter :: sigSB=5.6704d-8
  real(8) flux2T

  flux2T = ((1-albedo)*Q/sigSB/emiss)**0.25
end function flux2T


function interp1(x0,x,y0,y,xi,nz)
  ! linear interpolation
  ! x0<x(1)<x(2)<...<x(nz)
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: x(nz),y(nz),xi,x0,y0
  integer i
  real(8) interp1, yi
  
  yi = -9.
  do i=1,nz-1
     if (x(i)>x(i+1)) error stop 'incorrect direction'
     if (x(i)<xi) then
        yi = (y(i)*(x(i+1)-xi)+y(i+1)*(xi-x(i)))/(x(i+1)-x(i))
     endif
  enddo
  if (xi<x(1)) yi = (y0*x(1)-y(1)*x0)/(x(1)-x0)
  
  interp1 = yi
end function interp1


function heatcapacity(T)
  ! specific heat capacity of silicates
  implicit none
  real(8), intent(IN) :: T  ! [K]
  real(8) heatcapacity  ! J/(kg K)
  real(8) c 
  
  ! heat capacity from Ledlow et al. (1992), <350K
  !c = 0.1812 + 0.1191*(T/300.-1) + 0.0176*(T/300.-1)**2 + &
  !     0.2721*(T/300.-1)**3 + 0.1869*(T/300.-1)**4
  !c = c*1000*4.184  ! cal/(g K) -> J/(kg K)

  ! heat capacity from Winter & Saari (1969),  20K<T<500K
  c = -0.034*T**0.5 + 0.008*T - 0.0002*T**1.5
  heatcapacity = c*1000   ! J/(g K) -> J/(kg K)
end function heatcapacity


function vapordiffusivity(diam,porosity,T)
  ! diam = rms grain diameter
  implicit none
  real*8, parameter :: pi=3.1415926535897932
  real*8, parameter :: Ru = 8314.5
  real(8) vapordiffusivity
  real(8), intent(IN) :: diam,porosity,T
  real(8) vbar, r
  real(8), parameter :: tau = 2.  ! tortuosity

  r = diam/2.
  vbar = sqrt(8*Ru*T/(pi*18))
  ! for 0<=porosity<=0.5
  vapordiffusivity = pi/(8+pi)*porosity/(1-porosity)*vbar*r/tau
end function vapordiffusivity


function faintsun(t)
  implicit none
  real(8) faintsun
  real(8), intent(IN) :: t   ! time before present [years]
  ! Gough, D. O. (1981), Sol. Phys., 74, 21–34
  faintsun = 1./(1+0.4*abs(t)/4.57e9)
end function faintsun


