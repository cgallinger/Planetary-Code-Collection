!**********************************************************************
! Regolith thermal properties for the Moon
! Options include those used in the "Diviner standard thermal model"
! described in Hayne et al., JGR 122, 2371 (2017)
!**********************************************************************


function conductivity(rho, T)
  ! thermal conduction through solid
  implicit none
  real(8) conductivity
  real(8) Aam,Bam,Cam,Dam,Eam,Fam,Gam,Ham,Iam,kam
  real(8) A1,A2,B1,B2
  real(8), intent(IN) :: rho, T
  
  ! Thermal conductivity model from Martinez et al. (2021):
  
  ! Constants for temperature-dependent amorphous solid conductivity matrix found in Woods-Robinson et el. (2019)
  Aam = -2.03297e-1
  Bam = -11.472
  Cam = 22.5793
  Dam = -14.3084
  Eam = 3.41742
  Fam = 0.01101
  Gam = -2.80491e-5
  Ham = 3.35837e-8
  Iam = -1.40021e-11
  kam = Aam + Bam*(T**(-4)) + Cam*(T**(-3)) + Dam*(T**(-2)) + Eam*(T**(-1)) + Fam*(T) + Gam*(T**2) + Ham*(T**3) + Iam*(T**4)

  ! Derived constants for density-dependent functions
  A1 = (5.0821e-6)
  A2 = (-0.0051)
  B1 = 2.022e-13
  B2 = -1.953e-10

  conductivity = ((A1*rho + A2))*kam + (B1*rho + B2)*T**3

end function conductivity



pure function heatcapacity(T)
  ! specific heat capacity of silicates
  implicit none
  real(8), intent(IN) :: T  ! [K]
  real(8) c, heatcapacity   ! [J/(kg K)]

  ! Ledlow et al., ApJ 348, 640 (1992), <350K
  !c = 0.1812 + 0.1191*(T/300.-1) + 0.0176*(T/300.-1)**2 + &
  !     0.2721*(T/300.-1)**3 + 0.1869*(T/300.-1)**4
  !heatcapacity = c*1000*4.184  ! cal/(g K) -> J/(kg K)
  
  ! Winter & Saari, ApJ 156, 1135 (1969), 20K<T<500K
  !c = -0.034*T**0.5 + 0.008*T - 0.0002*T**1.5
  !heatcapacity = c*1000   ! J/(g K) -> J/(kg K)

  ! Hayne et al., JGR 122, 2371 (2017)
  c = 8.9093E-09*T**4 -1.2340E-05*T**3 +2.36160E-03*T**2 + 2.7431*T -3.6125
  heatcapacity = c

  ! Hemingway et al., Proc. 4th Lun. Sci. Conf. 3, 2481-2487 (1973), 90K<T<300K
  !c = —2.3173e-2 + 2.1270e-3*T + 1.5009e-5*T**2 — 7.3699e-8*T**3 + 9.6552e-11*T**4
  !heatcapacity = c

  ! Biele et al., Int. J. Thermophys. 43, 144 (2022), eq. 24
  !x = log(T)
  !c = exp((3.*x**3 -54.45*x**2 +306.8*x -376.6)/(x**2 -16.81*x +87.32))
  !heatcapacity = c
end function heatcapacity



function radconductivity(T, ell, emiss, porosity)
  ! radiative contribution to thermal conductivity
  implicit none
  real(8) radconductivity
  real(8), parameter :: sigSB = 5.6704d-8
  real(8), intent(IN) :: T, ell, emiss, porosity
  real(8) A
  
  A = 4 * emiss/(2-emiss) * sigSB * ell * (porosity/(1-porosity))**(1./3.)
  ! A = kc*chi/350**3

  ! Sakatani et al. (2017)
  ! zeta = 0.12
  ! A = 8*emiss/(2-emiss)*sigSB*zeta*( porosity/(1-porosity) )**(1./3.)*ell/2

  ! Gundlach & Blum (2013)
  ! A = 8*emiss*sigSB*1.34* porosity/(1-porosity) *ell/2
  
  radconductivity = A*T**3  ! also known as kr
end function radconductivity



function radconductivity1(T,chi,kc)
  ! radiative contribution to thermal conductivity
  implicit none
  real(8) radconductivity1
  real(8), intent(IN) :: T, chi, kc
  radconductivity1 = kc * ( 1 + chi*(T/350.)**3 )
end function radconductivity1



pure function twolayers(s,d,z,H)
  ! s ... surface value (z=0)
  ! d ... value at great depth (z=Inf)
  ! z ... depth
  ! H ... depth-scale
  implicit none
  real(8) twolayers
  real(8), intent(IN) :: s, d, z, H
  !H = 0.07 ! Hayne et al. (2017) average
  !H = 0.06 ! Hayne et al. (2017) Fig. A2
  twolayers = d - (d-s) * exp(-z/H)

  ! rho = (1-porosity)*rhosolid,  porosity = 1-rho/rhosolid
  ! e.g., rho_s = 1100; rho_d = 1800
  ! 1-1100/2500 = 0.560 ! s
  ! 1-1800/2500 = 0.280 ! d
  ! porosity(z) = 0.28 + (0.56-0.28)*exp(-z/H)
end function twolayers
