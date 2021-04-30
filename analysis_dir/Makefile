FC = gfortran

#FFLAGS = -O3 -g -fbacktrace
FFLAGS = -O3 -fbounds-check

PROGRAM = analysis

MODULES = module.f90
SOURCE = analysis.f90 readin.f90 calculate.f90

OBJECTS = $(SOURCE:.f90=.o)
MOD = $(MODULES:.f90=.o)


$(PROGRAM): $(OBJECTS) $(MOD) 
	$(FC) $(FFLAGS) $(MOD) $(OBJECTS) -o $(PROGRAM)


%.o : %.f90 $(MOD)
	$(FC) $(FFLAGS) -c $<
	
$(MOD): $(MODULES)
	$(FC) $(FFLAGS) -c $(MODULES)
	
clean :
	rm *.o *.mod $(PROGRAM)
