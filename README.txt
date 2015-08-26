Planetary Code Collection

by Norbert Schorghofer


This program collection contains
* a semi-implicit one-dimensional thermal model for planetary surfaces,
* an explicit subsurface vapor diffusion and deposition model,
* models for the equilibrium ice table on Mars,
* fast method for subsurface-atmosphere vapor exchange for Mars,
* a Monte-Carlo model for ballistic hops in the lunar exosphere, and
* a 3D model of shadowing, illumination, and scattering.



==Mars Subsurface Ice Model (M-SIM)==

*Mars Thermal Model*

Mars/mars_thermal1d.f: (main program)
Mars/flux.f
Mars/marsorbit.f
Common/conductionQ.f
Common/conductionT.f
Common/tridag.for
Common/grids.f
Common/psvco2.f
Mars/input.par
Documentation: Common/modeldescription.pdf (Part 1)


*Vapor Diffusion Model*

Mars/exper_thermal1d.f: (main program)
Mars/vapordiffusioni.f
Mars/adsorption.f
Mars/exper.par
Documentation: Common/modeldescription.pdf (Part 2)
Documentation: Schorghofer, N. & Aharonson, O. J. Geophys. Res. 110, E05003 (2005), Appendices


*Equilibrium Ice Table (on Mars)*

Mars/mars_mapi.f: (main program)
Mars/mars_mapt2.f: (main program)
Mars/mars_mapiqh2v.f90: (main program)
Mars/jsub.f
Mars/jsubv.f90
Mars/marsorbit.f
Common/conductionQ.f
Common/conductionT.f
Common/tridag.for
Common/grids.f
Common/julday.for
Common/psv.f
Common/psvco2.f
Mars/mapgrid.dat


*Fast Method for Subsurface Ice Dynamics*

Mars/stabgrow_fast.f90: (main program)
Mars/exper_fast.f90: (main program)
Mars/mars_fast.f90: (main program)
Mars/fast_modules.f90
Mars/fast_subs_univ.f90
Mars/fast_subs_exper.f90
Mars/fast_subs_mars.f90
Common/conductionQ.f
Common/conductionT.f
Common/tridag.for
Common/grids.f
Common/generalorbit.f
Common/psv.f
Common/psvco2.f
Common/derivs.f90
Mars/input_fast.par
Mars/lats.ph
Documentation: Schorghofer, N. Icarus 208, 598-607 (2010)


==Other Models for Planetary Surfaces==

*Basic Thermal Model for Asteroids*

Asteroid/asteroid_thermal.f90: (main program)
Asteroid/oneasteroid_thermal1d.f90
Asteroid/insolonly.f90 
Common/flux_noatm.f90
Common/conductionQ.f
Common/tridag.for
Common/grids.f
Common/generalorbit.f
Documentation: Common/modeldescription.pdf (Part 1)


*Monte-Carlo Model for Surface-bounded Exosphere*

Exosphere/moon4.f90: (main program)
Exosphere/montecarlo.f90
Exosphere/geogrid.f90
Exosphere/geogrid_D.f90
Common/subl_subs.f90
Common/gasdev.for
Common/ran2.for
Documentation: None


*Shadowing and illumination*

Craters/shadows.f90: (main program)
Craters/fieldofviews.f90: (main program)
Craters/cratersQ_snapshot.f90: (main program)
Craters/cratersQ_moon.f90: (main program)
Craters/topos.f90
Craters/crater_modules.f90
Craters/crater_common.f90
Craters/shadow_subs.f90
Craters/fieldofview_subs.f90
Craters/model_subs.f90
Common/hpsort.for
Craters/topo40.xyz
Documentation: None



ACKNOWLEDGMENTS:

2006: Troy Hudson discovered a grid-point offset in conductionT and conductionQ, which has been corrected.

2005: Thanks to Mischa Kreslavsky for providing correct formulas for energy balance on a slope.

* Many Thanks to Andy Vaught for developing an open-source Fortran 95 compiler (www.g95.org).  The early versions of this code were developed with this compiler.

2001: Samar Khatiwala invented an elegant implementation of the upper radiation boundary condition for the Crank-Nicolson method.

SUPPORT: This code development was supported by NASA, Caltech, and the University of Hawaii. Undoubtedly, some parts were written without support.


