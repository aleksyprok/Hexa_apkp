module grid
   use var_global
   implicit none

   contains

  subroutine grid_setup
    integer :: num_hex_cells,i,j,k,imax,jmax
    real    :: dmax,xoff,yoff
    integer, dimension(mpidir) :: dumcord
  integer :: namelen
  character*(MPI_MAX_PROCESSOR_NAME) :: procname


!
! Determine sub array division of global arrays.
! USE MPI routine to split up  grid - dims gives number of processes in
! each di rection.

  call MPI_DIMS_CREATE(mpisize,mpidir,dims,ierr)
  nproc = dims
  print*,mpisize,mpidir,dims,ierr,nproc

! USE MPI routine to define cartesian grid.

!   if (mpidir .eq. 1) then
!    periods(1) = .true.        ! periodic domain.
!   else
!    periods(1) = .false.
!    periods(2) = .false.
!   endif

  if (periodic .eq. 1) then
    periods(1)=.true.
    periods(2)=.true.
    print *, 'PERIODIC'
  else
    periods(1)=.false.
    periods(2)=.false.
    print *, 'CLOSED'
  endif


! sets of cartesian geometry with boundaries and places communicator in index comm.

  call MPI_CART_CREATE(MPI_COMM_WORLD,mpidir,nproc,periods,.TRUE.,comm,ierr)

  call MPI_COMM_RANK(comm,rank,ierr)

  CALL MPI_CART_COORDS(comm,rank,mpidir,coords,ierr)

  call MPI_Get_processor_name(procname,namelen,ierr)


  print *, 'rank',rank, procname(1:namelen),coords(1),coords(2)
  if (coords(1) .eq. 2) then
     print*,rank,coords,nproc(1)-1
  endif
  call MPI_BARRIER(comm,ierr)


!
! Determine rank of processes left/right and up/down of present rank
!

  CALL MPI_CART_SHIFT(comm,0,1,left,right,ierr)
  print*,'Left/right',left,rank,right
  if (mpidir .eq. 2) then
    CALL MPI_CART_SHIFT(comm,1,1,down,up,ierr)
    print*,'up down',down,rank,up
  endif


  dumcord = 0
  call MPI_CART_RANK(comm,dumcord,rankstart,ierr) !get rank of process at (0,0)


!
!  Read param1 file and compute spacings in 3d array:
!
    if (rank .eq.rankstart) then
      filename=dir(1:length(dir,20))//'/param1'
                                           !     print *,'Read grid parameters: ',filename
      open(unit=3,file=filename,form='formatted',status='old')
      read (3,*) num_hex_cells
      read (3,*) i,j,k
      read(3,*) length_cm
      read(3,*) time_s
      close(3)
    endif
    CALL MPI_BCAST(num_hex_cells,1,MPI_INTEGER,rankstart,comm,ierr)
    CALL MPI_BCAST(i,1,MPI_INTEGER,rankstart,comm,ierr)
    CALL MPI_BCAST(j,1,MPI_INTEGER,rankstart,comm,ierr)
    CALL MPI_BCAST(k,1,MPI_INTEGER,rankstart,comm,ierr)
    CALL MPI_BCAST(length_cm,1,MPI_real,rankstart,comm,ierr)
    CALL MPI_BCAST(time_s,1,MPI_real,rankstart,comm,ierr)
                                                      !

    call MPI_BARRIER(comm,ierr)

    nxglobal=i
    nyglobal=j
    nzglobal=k

    if( (i.ne.nxglobal) .or. (j.ne.nyglobal) .or. (k.ne.nzglobal) ) then
       print*,"arrays do not match theoretical limit : correct",rank
       call MPI_FINALIZE(ierr)
       stop
    endif

   nz = nzglobal                            ! no mpi in z direction
 if (mpidir .eq. 1) then
   nx = (nxglobal)/nproc(1)             ! no of cells associated with each process.
   ny = nyglobal                          ! no mpi in y direction
                                         !   print*,nx,ny,nz

   if ( (nx*nproc(1)) .ne. nxglobal) then
     print*,'Unable to subdivide equally. Fix grid'
     CALL MPI_FINALIZE(ierr)
     stop
   endif

 else
   nx = (nxglobal)/nproc(1)            ! no of cells in x associated with each process.
   ny = (nyglobal)/nproc(2)            ! no of cells in y associated with each process.
                                       !   print*,nx,ny,nz

   if ( ((nx*nproc(1)) .ne. nxglobal) .or.  ((ny*nproc(2)) .ne. nyglobal) ) then
       if ( (nx*nproc(1)) .ne. nxglobal) then
          print*,'Unable to subdevide equally in x. Fix grid'
       endif
       if (  (ny*nproc(2)) .ne. nyglobal ) then
          print*,'Unable to subdivide equally in y. Fix grid'
       endif
       CALL MPI_FINALIZE(ierr)
       stop
   endif

 endif


!
!  Compute spacings in 3D array:
!

    xoff=0.5
    yoff=0.5
    imax=num_hex_cells
    jmax=nint(imax*sqrt(3.0))
    dmax=3.0*imax-1.0+2.*xoff
    delx=dmax/nxglobal
    dely=dmax/nyglobal
    delz=delx
                                     !    print *,'delx=',delx,', dely=',dely,', delz=',delz


  end subroutine grid_setup

!-----------------------------------------------------------------------
! Enter paramter values for runs.
!-----------------------------------------------------------------------

subroutine setup_param


    if (rank .eq. 0) then
      write (filename,fmt='(a,"/",a,"_setup")') &
         dir(1:length(dir,20)),root(1:length(root,10))
      print *,'Read model parameters: ',filename
      open(unit=3,file=filename,form='formatted',status='old')
      vsetup=get_value(3,'vsetup')    ! setup file version
      nmajor=get_value(3,'nmajor')
      nstrt  =get_value(3,'nstrt')
      nend  =get_value(3,'nend')
      etaia  =get_value(3,'etaia')
      eta4a  =get_value(3,'eta4a')
      periodic =get_value(3,'periodic')
      open=get_value(3,'open')

      close(3)

      print *,'nmajor=',nmajor
      print *,'nstrt=',nstrt
      print *,'nend =',nend
      print *,'etaia =',etaia
      print *,'eta4a =',eta4a
      print *,'periodic=',periodic
      print *, 'open =',open

   endif

  CALL MPI_BCAST(vsetup,1,MPI_INTEGER,0,comm,ierr)
  CALL MPI_BCAST(nmajor,1,MPI_INTEGER,0,comm,ierr)
  CALL MPI_BCAST(nminor,1,MPI_INTEGER,0,comm,ierr)

  CALL MPI_BCAST(nstrt,1,MPI_INTEGER,0,comm,ierr)
  CALL MPI_BCAST(nend,1,MPI_INTEGER,0,comm,ierr)
  CALL MPI_BCAST(etaia,1,MPI_REAL,0,comm,ierr)
  CALL MPI_BCAST(eta4a,1,MPI_REAL,0,comm,ierr)
  CALL MPI_BCAST(periodic,1,MPI_INTEGER,0,comm,ierr)
  CALL MPI_BCAST(open,1,MPI_INTEGER,0,comm,ierr)

end subroutine setup_param

  subroutine arrayaloc

    allocate(aax(0:nx+1,0:ny+2,nz+1),aay(0:nx+2,0:ny+1,nz+1),aaz(0:nx+2,0:ny+2,nz))
    allocate(aax0(0:nx+1,0:ny+2,ntmax),aay0(0:nx+2,0:ny+1,ntmax))
    allocate(daax(0:nx+1,0:ny+2),daay(0:nx+2,0:ny+1))

    allocate(bbx(nx+1,ny+2,nz+2),bby(nx+2,ny+1,nz+2),bbz(nx+2,ny+2,nz+1))
    allocate(ccx(nx+2,ny+1,nz+1),ccy(nx+1,ny+2,nz+1),ccz(nx+1,ny+1,nz+2))
    allocate(div(nx+1,ny+1,nz+1),bx(nx+1,ny+1,nz+1),by(nx+1,ny+1,nz+1),bz(nx+1,ny+1,nz+1))
    allocate(bb(nx+1,ny+1,nz+1),bbm(nx+1,ny+1,nz+1),ch(nx+1,ny+1,nz+1))
    allocate(cx(nx+1,ny+1,nz+1),cy(nx+1,ny+1,nz+1),cz(nx+1,ny+1,nz+1))
    allocate(vx(nx+1,ny+1,nz+1),vy(nx+1,ny+1,nz+1),vz(nx+1,ny+1,nz+1))

    if (open .eq. 1) then !allocate flux correction array
      allocate(ay_corr(0:nx+2,0:ny+1,nz+1))
    endif

  end subroutine arrayaloc



  subroutine arraydealoc

    deallocate(aax,aay,aaz)
    deallocate(aax0,aay0)
    deallocate(daax,daay)
    deallocate(bbx,bby,bbz)
    deallocate(ccx,ccy,ccz)
    deallocate(div,bx,by,bz)
    deallocate(bb,bbm,ch)
    deallocate(cx,cy,cz)
    deallocate(vx,vy,vz)

  end subroutine arraydealoc


endmodule grid
