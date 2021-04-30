program corona
!creates a linear force-free (or potential) initial condition coronal field
!for use in Hexa. This code reads in the 'params.dat' file (which contains the
!gridsize, working directory and alpha value) and 'mag_data' (containing the
!time series of cleaned magnetograms) produced by 'data_prep.pro' and produces
! the '_00000p' file.

use avalues_mod
use dunk_mod
use sysdat_mod
use cor_mag_mod



      print*,'#########################################'
      print*,'######### Coronal field generator #######'
      print*,'#########################################'
      print*,''


      CALL sysdat !Read in 'params.dat' and setup gridsize etc.

      CALL array_alloc() !allocate arrays

      CALL readin ! readin 'mag_data' file

      allocate(ax(2*na+1,2*na+1),ay(2*na+1,2*na+1))




      dunk_count=1

      call poislin() ! calculate Ax and Ay from Bz of first magnetogram at photosphere

      call array_dealloc() !deallocate evolve arrays
      deallocate(cor_mag)

       print*,''

      Print*,'Generating LFF field.'
      if (2*na .eq. 256) then
        print*,'This should take around 1 minute for the chosen gridsize'
      else if (2*na .eq. 512) then
        print*,'This should take around 10 minutes for the chosen gridsize'
      else
        print*,'This should take a few seconds for the chosen gridsize'
      endif
      call f3cray() !generates 3d linear force-free field

      print*,''
      call atob() !takes A, averages it onto cell ribs and writes to '_00000p'

      end
