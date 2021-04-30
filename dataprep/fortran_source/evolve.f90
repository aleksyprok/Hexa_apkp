program surface_evol
!This code produces the evolving lower boundary condition used by Hexa.
!The code reads in the 'params.dat' file (which contains the gridsize, working
!directory and number of magnetogram frames) and 'mag_data' (containing the time
!series of cleaned magnetograms) produced by 'data_prep.pro'.
!This code outputs the '_evolve' file

use avalues_mod
use dunk_mod
use sysdat_mod
use fnames


print*,'#########################################'
print*,'###### Surface evolution generator ######'
print*,'#########################################'
print*,''



      CALL sysdat !get gridsize, number of frames and working directory

      CALL array_alloc() !allocate arrays

      CALL readin ! readin 'mag_data' file


      allocate(ax(2*na+1,2*na+1),ay(2*na+1,2*na+1))


      open (unit=11,file=trim(dir)//'/hexa_files/'//trim(dir)//'_evolve',form='unformatted')

      j = 1

      write(11) 1


      dunk_count=0
      do i=1,nt !loop over frames
         dunk_count = dunk_count+1
         call poislin() ! calculate Ax and Ay from Bz on the base
         call writeout_evolve !write Ax and Ay to the '_evolve' file
         write(*,"('Frame ',i3,' of ',i3,' processed')")dunk_count,nt
      enddo
      print*,''
      print*,'File '//trim(dir)//'/hexa_files/'//trim(dir)//'_evolve written to file.'
      call array_dealloc

      close(11)


      end
