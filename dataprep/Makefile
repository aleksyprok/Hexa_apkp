compiler = gfortran -O3 -ffast-math

modules = fortran_source/modules.f90
evolve_subs = fortran_source/poislin.f90 fortran_source/fluxfinder.f90
corona_subs = fortran_source/f3cray.f90

all: 
	$(compiler) $(modules) $(evolve_subs) fortran_source/evolve.f90 -o evolve
	$(compiler) $(modules) $(evolve_subs) $(corona_subs) fortran_source/corona.f90 -o corona
	rm *.mod

evolve: $(modules) $(evolve_subs) fortran_source/evolve.f90
	$(compiler) $(modules) $(evolve_subs) fortran_source/evolve.f90 -o evolve
	rm *.mod

corona: $(modules) $(evolve_subs) $(corona_subs) fortran_source/corona.f90
	$(compiler)  $(modules) $(evolve_subs) $(corona_subs) fortran_source/corona.f90 -o corona
	rm *.mod
	
