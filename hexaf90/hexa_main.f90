program hexa
  use var_global
  use grid
  use io
  use cal
  implicit none

  integer :: nt,it,i,j,k,kmax,ios,n,nr,unlink
  integer :: num_steps,nma,nma1,nma2,nmi,itmax
  real    :: dmax,dy
  character*50 :: dummy
!
! set filename and directory
!
  dir='run1'
  root='run1'

  timestep_s=60. !timestep in seconds. This is used to define the diffusion coefficients.

  CALL MPI_INIT(ierr) ! initiate MPI

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,ierr) ! get # of processes

!
! set run parameters
!

  call setup_param !reads in _setup file

  call grid_setup !sets up computational grid and MPI (and reads in param1)



  basedt=timestep_s/time_s !timestep used to define frictional coefficient
  dt=minval( (/ basedt,0.1 /) ) !set dt to be basedt initially or 0.1 of the magnetogram cadence (whichever is smaller)
  ntmax=nmajor

  if (open .eq. 1) then
    allocate(corr(ntmax-1)) !allocate flux imbalance correction if open boundary condition
    allocate(fimb(ntmax))
  endif



  call arrayaloc !allocates arrays

  call MPI_Barrier(MPI_comm_world,ierr)




  etad= 0.05*delx**2/basedt

  etai=etaia/length_cm**2*time_s*1.e10
  eta4=eta4a/length_cm**4*time_s*1.e10

  if (rank .eq. 0) then
     print*, 'Diffusion of div A constant=',etad*length_cm**2/time_s/(1.e5**2),'km^2/s'
     print*, 'Ohmic diffusion constant=',etai*length_cm**2/time_s/(1.e5**2),'km^2/s'
     print*, 'hyperdiffusion constant=',eta4*length_cm**4/time_s/(1.e5**4),'km^4/s'
  endif

!
!  Read 3d vector potentials from file:
!

  nt=nstrt
  time=nstrt
  call readdata(nt,nr,'3d')    ! read in full 3d data

!
!  Evolve 3D model with time-dependent boundary conditions:
!
  call readdata(ntmax,nr,'2d') !read in lower boundary condition


  nma1=nt     ! start frame
  nma2=min(nend,nmajor-1) !stop frame
  if (rank .eq. rankstart) then
     open (unit=50,file='diagnostic',status='unknown')
  endif


!
!   Correction to Ay ghost cells if open top boundary and periodic side boundaries are used
   if (periodic .eq. 1 .and. open .eq. 1) then
	    if (coords(1) .eq. 0) then !left most processors
		aay(0,:,:)=aay(0,:,:)-fimb(nma1+1)*real(nxglobal)*delx
	    endif
	    if (coords(1) .eq. nproc(1)-1) then !right most processors
	      aay(nx+1,:,:)=aay(nx+1,:,:) +fimb(nma1+1)*real(nxglobal)*delx
	      aay(nx+2,:,:)=aay(nx+2,:,:) +fimb(nma1+1)*real(nxglobal)*delx
	    endif
   endif
!


  do nma=nma1,nma2 -1
	 !reminder: hexa files start from 0 and fortran arrays start from 1
	  n=nma+1 !hence why n=nma+1 ^

	  time=real(nma)

	  if (open .eq. 1) then
	     imb=corr(n) !flux imbalance correction for current frame
	     totimb=fimb(n)
	  endif



	  daax(0:nx+1,0:ny+2)=(aax0(0:nx+1,0:ny+2,n+1)-aax0(0:nx+1,0:ny+2,n)) !dAx/dt
	  daay(0:nx+2,0:ny+1)=(aay0(0:nx+2,0:ny+1,n+1)-aay0(0:nx+2,0:ny+1,n)) !dAy/dt


	  call step!(mode,nrep,'y')

	  nt=nt+1

	  call writedata(n,nr,'3d')

	  if (rank .eq. rankstart) then
		write (filename,fmt='(a,"/",a,"_log")') &
		    dir(1:length(dir,30)),root(1:length(root,10))
		print *,'Update log file: ',filename
		open(unit=3,file=filename,form='formatted',status='unknown',position='append')
		write (3,99) mode,nstrt,nend,  &
		    etaia,eta4a
		close(3)
		print *,'OK'
	  endif

  enddo

  call MPI_BARRIER(comm,ierr)


  call arraydealoc
  call MPI_Finalize(ierr)

  99        format("mode=",i1,", nstrt=",i5,", nend=",i5,  &
                           ", etai=",f7.5,", eta4=",f7.5)


end program hexa
