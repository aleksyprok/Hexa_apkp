compiler = mpif90 -O3 -ffast-math

#compiler = mpiifort -O3 -traceback -gen-interfaces -warn interfaces -check
#compiler=mpiifort -O3

modules = modules.f90
evolve_subs = poislin.f90 fluxfinder.f90
corona_subs = f3cray.f90

# all: 
# 	$(compiler) $(modules) $(evolve_subs) evolve.f90 -o evolve
# 	$(compiler) $(modules) $(evolve_subs) $(corona_subs) corona.f90 -o corona
# 	rm *.mod
# 
# evolve: $(modules) $(evolve_subs) evolve.f90
# 	$(compiler) $(modules) $(evolve_subs) evolve.f90 -o evolve
# 	rm *.mod

corona: $(modules) $(evolve_subs) $(corona_subs) corona.f90
	$(compiler)  $(modules) $(evolve_subs) $(corona_subs) corona.f90 -o corona
	rm *.mod
	
