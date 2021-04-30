FFLAGS =#F90 = mpiifort #-O3 #-xHost 
F90 = mpif90 -O3 --free-line-length-400 #-fbounds-check
#F90 = /share/apps/bin/mpif90 -O3
#F90 = mpif90 -O3 --free-line-length-400 -frecord-marker=4        
#F90 = mpif90 -g -trapuv -Wunitialized -O2
#F90 =  mpif90 -f90=pgf90 -Ktrap=fp -Mchkfpstk -O2

FFILES1 = hexa_var.o hexa_grid.o hexa_io.o hexa_cal.o hexa_main.o

hexa:$(FFILES1)
	$(F90) $(FFLAGS) $(FFILES1) -o hexa

hexa_var.o:hexa_var.f90
	$(F90) -c hexa_var.f90 $(FFLAGS1)

hexa_grid.o:hexa_grid.f90 hexa_var.o
	$(F90) -c hexa_grid.f90 $(FFLAGS1)

hexa_io.o:hexa_io.f90 hexa_var.o
	$(F90) -c hexa_io.f90 $(FFLAGS1)

hexa_cal.o:hexa_cal.f90 hexa_var.o
	$(F90) -c hexa_cal.f90 $(FFLAGS1)

hexa_main.o:hexa_main.f90 hexa_var.o hexa_grid.o hexa_io.o hexa_cal.o
	$(F90) -c  hexa_main.f90 $(FFLAGS1)

clean:
	rm *% *.o *.mod core *~ *gtest* *out.dat* hexa
