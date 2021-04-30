program corona
use MPI
use avalues_mod
use dunk_mod
use sysdat_mod
use cor_mag_mod
use cdata
use fnames
    
      integer :: comm, rank, size
      character (len=60) :: fmt
    
      !Initialise MPI, get rank and size
      call MPI_Init(ierror)      
      comm=MPI_COMM_WORLD      
      call MPI_Comm_size(comm,size,ierror)     
      call MPI_comm_rank(comm,rank,ierror)
      
      
      
      if (rank .eq. 0) then !print welcome message and read in params.dat
        print*,'#########################################'
        print*,'###### MPI Potential field generator ####'
        print*,'#########################################'
        print*,''       
        CALL sysdat !get gridsize
        if (alpha .ne. 0.) then
          print*, "###########################"
          print*, "Warning, Alpha is not zero!"
          print*, "Setting Alpha to zero now"
          print*, "###########################"
          alpha=0.d0 !force alpha to be zero
        endif
        print*,'' 
        fmt='("Number of processors=", i3)'
        write(*,fmt),size
        print*,''        
      endif
      
      call transmit !broadcast parameters to other processes
      
      CALL array_alloc() !allocate arrays
      
      if (rank .eq. 0) then
	  CALL readin ! readin 'mag_data' file 
      endif 
      
      call MPI_Bcast(cor_mag,2*na*2*na*nt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
      
      allocate(ax(2*na+1,2*na+1),ay(2*na+1,2*na+1))
      
      
      
      !count the number of frames to process per processor
      nperproc=ceiling(real(nt)/real(size))
      
      
      
      itop=(rank+1)*nperproc
      if (rank .eq. size-1) itop=nt
      ibottom=rank*nperproc+1
      
      call MPI_barrier(comm,ierror)
      fmt='("Processor ",I2," processes frames ",i3," to ",i3)'
      !print*,'Processor ',rank+1,'processes frames',ibottom,'to',itop
      write(*,fmt) rank+1,ibottom,itop
     
      
      call MPI_barrier(comm,ierror)
      
      
      do i=ibottom,itop
	  dunk_count=i
	  
	  call poislin() ! calculate Ax and Ay from Bz on base        
	  call f3cray() !field solver

	  call atob() !takes a, and writes to poten_00000p
	  
      enddo
      
      Print*,'Processor',rank+1,' finished!'
      call MPI_Finalize(ierror)
      end

