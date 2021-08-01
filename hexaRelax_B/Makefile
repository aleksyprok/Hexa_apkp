F90 = mpif90 --free-line-length-400

FFILES1 = relax_var.o relax_grid.o relax_io.o relax_cal.o relax_main.o
# FFLAGS1 = -O0 -g -Wall -Wextra -pedantic -fbounds-check \
#           -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter \
# 			    -Wrealloc-lhs-all -fno-realloc-lhs
FFLAGS1 = -O3 -fallow-argument-mismatch

relax:$(FFILES1)
	$(F90) $(FFLAGS1) $(FFILES1) -o relax

relax_var.o:relax_var.f90
	$(F90) -c relax_var.f90 $(FFLAGS1)

relax_grid.o:relax_grid.f90 relax_var.o
	$(F90) -c relax_grid.f90 $(FFLAGS1)

relax_io.o:relax_io.f90 relax_var.o
	$(F90) -c relax_io.f90 $(FFLAGS1)

relax_cal.o:relax_cal.f90 relax_var.o
	$(F90) -c relax_cal.f90 $(FFLAGS1)

relax_main.o:relax_main.f90 relax_var.o relax_grid.o relax_io.o  relax_cal.o
	$(F90) -c  relax_main.f90 $(FFLAGS1)

clean:
	rm -rf *.o *.mod relax Figures/* run1/*
