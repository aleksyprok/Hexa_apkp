module io
   use var_global
   implicit none

   contains

  subroutine readdata(ndum,nr,tp)

    character*2 :: tp
    integer, intent(in) ::ndum, nr
    integer :: i,j,k,n
    integer, dimension(mpidir) :: dumcord
    real, dimension(:,:,:), allocatable :: aax_global, aay_global, aaz_global
    real, dimension(:,:,:), allocatable :: aax0_global, aay0_global

!    print*,rank,coords,rankstart

    if (tp .eq. '3d') then

                                                               !          aax = 0.0
                                                               !          aay = 0.0
                                                               !          aaz = 0.0
      if (rank .eq. rankstart) then

         allocate(aax_global(nxglobal,nyglobal+1,nzglobal+1))
         allocate(aay_global(nxglobal+1,nyglobal,nzglobal+1))
         allocate(aaz_global(nxglobal+1,nyglobal+1,nzglobal))

         write (filename,fmt='(a,"/",a,"_",i5.5,"p")')  &     ! *** Mac
              dir(1:length(dir,30)),root(1:length(root,10)),ndum
         print *,'Reading 3D model from '//filename

         open(unit=1,file=filename,form='unformatted',status='old')
         read(1)opt
         print *,'File type: opt=',opt
         if(.not.(opt.eq.1)) stop 'Invalid option'
         read(1)(((aax_global(i,j,k),i=1,nxglobal  ),j=1,nyglobal+1),k=1,nzglobal+1)
         read(1)(((aay_global(i,j,k),i=1,nxglobal+1),j=1,nyglobal  ),k=1,nzglobal+1)
         read(1)(((aaz_global(i,j,k),i=1,nxglobal+1),j=1,nyglobal+1),k=1,nzglobal  )
         close(1)

      endif
                                                              ! call MPI_BARRIER(comm,ierr)
!
! Transfer out to all processes.
!
      if (rank .eq. rankstart ) then

          do j=0,nproc(2)-1

            do i =0,nproc(1)-1

              dumcord(1) = i
              dumcord(2) = j
              call MPI_CART_RANK(comm,dumcord,nextrank,ierr)

              if ( ( i .eq. 0) .and. (j .eq. 0)) then
                  nextrank = mpi_proc_null
                  aax(1:nx,1:ny+1,1:nz+1) = aax_global(1:nx,1:ny+1,1:nz+1)
                  aay(1:nx+1,1:ny,1:nz+1) = aay_global(1:nx+1,1:ny,1:nz+1)
                  aaz(1:nx+1,1:ny+1,1:nz)   = aaz_global(1:nx+1,1:ny+1,1:nz)
              endif

	      if (j .ne. nproc(2)-1) then
                if ( ( i .ne. nproc(1)-1 ) ) then
                  CALL MPI_SEND(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aaz_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz),nx*ny*nz,MPI_REAL,nextrank,tag,comm,ierr)
                else
                  CALL MPI_SEND(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aaz_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:nz),(nx+1)*ny*nz,MPI_REAL,nextrank,tag,comm,ierr)
                endif
              else
                if ( (i .ne. nproc(1)-1 ) ) then
                  CALL MPI_SEND(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aaz_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:nz),nx*(ny+1)*nz,MPI_REAL,nextrank,tag,comm,ierr)
                else
                  CALL MPI_SEND(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aaz_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny+1,1:nz),(nx+1)*(ny+1)*nz,MPI_REAL,nextrank,tag,comm,ierr)
                endif
              endif

           enddo
         enddo
      else
          if (coords(1) .ne. nproc(1)-1) then
             if ( coords(2) .ne. nproc(2)-1) then
                 CALL MPI_RECV(aax(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aaz(1:nx,1:ny,1:nz),nx*ny*nz,MPI_REAL,rankstart,tag,comm,stat,ierr)
             else
                 CALL MPI_RECV(aax(1:nx,1:ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aaz(1:nx,1:ny+1,1:nz),nx*(ny+1)*nz,MPI_REAL,rankstart,tag,comm,stat,ierr)
             endif
          else
             if (coords(2) .ne. nproc(2)-1)  then
                 CALL MPI_RECV(aax(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay(1:nx+1,1:ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aaz(1:nx+1,1:ny,1:nz),(nx+1)*ny*nz,MPI_REAL,rankstart,tag,comm,stat,ierr)
             else
                 CALL MPI_RECV(aax(1:nx,1:ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay(1:nx+1,1:ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aaz(1:nx+1,1:ny+1,1:nz),(nx+1)*(ny+1)*nz,MPI_REAL,rankstart,tag,comm,stat,ierr)
             endif
          endif

      endif

CALL MPI_Barrier(comm,ierr) !Gordon: Ensure data is transferred before deallocation.

      if (rank .eq. rankstart) then
         deallocate(aax_global,aay_global,aaz_global)
      endif

!   do x-direction transfer
      if (coords(2) .ne. nproc(2)-1) then
          call MPI_SENDRECV(aax(nx,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,right,tag,aax(0,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aax(1,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,left,tag,aax(nx+1,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,right,tag,comm,stat,ierr)

          call MPI_SENDRECV(aay(nx,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,right,tag,aay(0,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aay(1:2,1:ny,1:nz+1),2*ny*(nz+1),MPI_REAL,left,tag,aay(nx+1:nx+2,1:ny,1:nz+1),2*ny*(nz+1),MPI_REAL,right,tag,comm,stat,ierr)

          call MPI_SENDRECV(aaz(nx,1:ny,1:nz),ny*nz,MPI_REAL,right,tag,aaz(0,1:ny,1:nz),ny*nz,MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aaz(1:2,1:ny,1:nz),2*ny*nz,MPI_REAL,left,tag,aaz(nx+1:nx+2,1:ny,1:nz),2*ny*nz,MPI_REAL,right,tag,comm,stat,ierr)
       else
          call MPI_SENDRECV(aax(nx,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,right,tag,aax(0,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aax(1,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,left,tag,aax(nx+1,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,right,tag,comm,stat,ierr)

          call MPI_SENDRECV(aay(nx,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,right,tag,aay(0,1:ny,1:nz+1),ny*(nz+1),MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aay(1:2,1:ny,1:nz+1),2*ny*(nz+1),MPI_REAL,left,tag,aay(nx+1:nx+2,1:ny,1:nz+1),2*ny*(nz+1),MPI_REAL,right,tag,comm,stat,ierr)

          call MPI_SENDRECV(aaz(nx,1:ny+1,1:nz),(ny+1)*nz,MPI_REAL,right,tag,aaz(0,1:ny+1,1:nz),(ny+1)*nz,MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aaz(1:2,1:ny+1,1:nz),2*(ny+1)*nz,MPI_REAL,left,tag,aaz(nx+1:nx+2,1:ny+1,1:nz),2*(ny+1)*nz,MPI_REAL,right,tag,comm,stat,ierr)
       endif

!   do vertical transfer.



       call MPI_SENDRECV(aax(0:nx+1,ny,1:nz+1),(nx+2)*(nz+1),MPI_REAL,up,tag,aax(0:nx+1,0,1:nz+1),(nx+2)*(nz+1),MPI_REAL,down,tag,comm,stat,ierr)
       call MPI_SENDRECV(aax(0:nx+1,1:2,1:nz+1),(nx+2)*2*(nz+1),MPI_REAL,down,tag,aax(0:nx+1,ny+1:ny+2,1:nz+1),(nx+2)*2*(nz+1),MPI_REAL,up,tag,comm,stat,ierr)

       call MPI_SENDRECV(aay(0:nx+2,ny,1:nz+1),(nx+3)*(nz+1),MPI_REAL,up,tag,aay(0:nx+2,0,1:nz+1),(nx+3)*(nz+1),MPI_REAL,down,tag,comm,stat,ierr)
       call MPI_SENDRECV(aay(0:nx+2,1,1:nz+1),(nx+3)*(nz+1),MPI_REAL,down,tag,aay(0:nx+2,ny+1,1:nz+1),(nx+3)*(nz+1),MPI_REAL,up,tag,comm,stat,ierr)

       call MPI_SENDRECV(aaz(0:nx+2,ny,1:nz),(nx+3)*nz,MPI_REAL,up,tag,aaz(0:nx+2,0,1:nz),(nx+3)*nz,MPI_REAL,down,tag,comm,stat,ierr)
       call MPI_SENDRECV(aaz(0:nx+2,1:2,1:nz),2*(nx+3)*nz,MPI_REAL,down,tag,aaz(0:nx+2,ny+1:ny+2,1:nz),2*(nx+3)*nz,MPI_REAL,up,tag,comm,stat,ierr)




       if (up .eq. MPI_PROC_NULL ) then
            aax(0:nx+1,ny+2,1:nz+1) = 0.0
            aay(0:nx+2,ny+1,1:nz+1) = 0.0
            aaz(0:nx+2,ny+2,1:nz) = 0.0
       endif

       if (down .eq. MPI_PROC_NULL) then
            aax(0:nx+1,0,1:nz+1) = 0.0
            aay(0:nx+2,0,1:nz+1) = 0.0
            aaz(0:nx+2,0,1:nz) = 0.0
       endif

       if (right .eq. MPI_PROC_NULL ) then
            aax(nx+1,0:ny+2,1:nz+1) = 0.0
            aay(nx+2,0:ny+1,1:nz+1) = 0.0
            aaz(nx+2,0:ny+2,1:nz) = 0.0
       endif

       if (left .eq. MPI_PROC_NULL) then
            aax(0,0:ny+2,1:nz+1) = 0.0
            aay(0,0:ny+1,1:nz+1) = 0.0
            aaz(0,0:ny+2,1:nz) = 0.0
       endif




!       print*,maxval(aax),minval(aax),maxloc(aax),minloc(aax)
!       print*,maxval(aay),minval(aay),maxloc(aay),minloc(aay)
!       print*,maxval(aaz),minval(aaz),maxloc(aaz),minloc(aaz)
endif

if (tp .eq.'2d') then
      aax0 = 0.0
      aay0 = 0.0

      if (rank .eq. rankstart) then

         allocate(aax0_global(nxglobal,nyglobal+1,ntmax))
         allocate(aay0_global(nxglobal+1,nyglobal,ntmax))

         write (filename,fmt='(a,"/",a,"_evolve")') &
              dir(1:length(dir,20)),root(1:length(root,10))
         print *,'Read boundary conditions: ',filename
         open(unit=1,file=filename,form='unformatted',status='old')
         read (1) opt
         print *,'File type: opt=',opt
         if(.not.(opt.eq.1)) stop 'Invalid option'
         do n=1,ndum
           print*,n
           read (1) ((aax0_global(i,j,n),i=1,nxglobal  ),j=1,nyglobal+1)
           read (1) ((aay0_global(i,j,n),i=1,nxglobal+1),j=1,nyglobal )
         enddo
        close(1)
        print *,'OK'


        !determine flux balance correction required per frame (if using open boundary conditions)
        if (open .eq. 1) then

           do i=1,ntmax ! flux = \oint A.dl
              fimb(i) = sum(aax0_global(:,1,i))*delx+sum(aay0_global(nxglobal+1,:,i))*dely &
                        -sum(aax0_global(:,nyglobal+1,i))*delx-sum(aay0_global(1,:,i))*dely
           enddo
            !fimb per pixel = flux / area / nx/ny
           fimb=(fimb/real(delx)/real(dely))/real(nxglobal)/real(nyglobal)

           !correction to be applied to each frame
           corr=fimb(2:ntmax)-fimb(1:ntmax-1)


        endif
      endif


      if (open .eq. 1) then !broadcast flux balance correction to all processes
        CALL MPI_BCAST(corr,ntmax-1,MPI_REAL,rankstart,comm,ierr)
        CALL MPI_BCAST(fimb,ntmax,MPI_REAL,rankstart,comm,ierr)
      endif
!
! Transfer out to all processes.
!
      if (rank .eq. rankstart ) then

          do j=0,nproc(2)-1

            do i =0,nproc(1)-1

              dumcord(1) = i
              dumcord(2) = j
              call MPI_CART_RANK(comm,dumcord,nextrank,ierr)

              if ( ( i .eq. 0) .and. (j .eq. 0)) then
                  nextrank = mpi_proc_null
                  aax0(1:nx,1:ny+1,1:ndum) = aax0_global(1:nx,1:ny+1,1:ndum)
                  aay0(1:nx+1,1:ny,1:ndum) = aay0_global(1:nx+1,1:ny,1:ndum)


              endif

	      if (j .ne. nproc(2)-1) then
                if ( ( i .ne. nproc(1)-1 ) ) then
                  CALL MPI_SEND(aax0_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:ndum),nx*ny*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay0_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:ndum),nx*ny*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                else
                  CALL MPI_SEND(aax0_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:ndum),nx*ny*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay0_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:ndum),(nx+1)*ny*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                endif
              else
                if ( (i .ne. nproc(1)-1 ) ) then
                  CALL MPI_SEND(aax0_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:ndum),nx*(ny+1)*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay0_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:ndum),nx*ny*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                else
                  CALL MPI_SEND(aax0_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:ndum),nx*(ny+1)*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                  CALL MPI_SEND(aay0_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:ndum),(nx+1)*ny*ndum,MPI_REAL,nextrank,tag,comm,ierr)
                endif
              endif

           enddo
         enddo
   else
          if (coords(1) .ne. nproc(1)-1) then
             if ( coords(2) .ne. nproc(2)-1) then
                 CALL MPI_RECV(aax0(1:nx,1:ny,1:ndum),nx*ny*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay0(1:nx,1:ny,1:ndum),nx*ny*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
             else
                 CALL MPI_RECV(aax0(1:nx,1:ny+1,1:ndum),nx*(ny+1)*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay0(1:nx,1:ny,1:ndum),nx*ny*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
             endif
          else
             if (coords(2) .ne. nproc(2)-1)  then
                 CALL MPI_RECV(aax0(1:nx,1:ny,1:ndum),nx*ny*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay0(1:nx+1,1:ny,1:ndum),(nx+1)*ny*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
             else
                 CALL MPI_RECV(aax0(1:nx,1:ny+1,1:ndum),nx*(ny+1)*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
                 CALL MPI_RECV(aay0(1:nx+1,1:ny,1:ndum),(nx+1)*ny*ndum,MPI_REAL,rankstart,tag,comm,stat,ierr)
             endif
          endif

      endif

call MPI_Barrier(comm,ierr) !Gordon - ensures all sending and recieving is done before deallocation

      if (rank .eq. rankstart) then
         print *, 'Barrier done. Deallocating global arrays.'
         deallocate(aax0_global,aay0_global)
      endif

         !   do x-direction transfer
      if (coords(2) .ne. nproc(2)-1) then
          call MPI_SENDRECV(aax0(nx,1:ny,1:ndum),ny*ndum,MPI_REAL,right,tag,aax0(0,1:ny,1:ndum),ny*ndum,MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aax0(1,1:ny,1:ndum),ny*ndum,MPI_REAL,left,tag,aax0(nx+1,1:ny,1:ndum),ny*ndum,MPI_REAL,right,tag,comm,stat,ierr)

          call MPI_SENDRECV(aay0(nx,1:ny,1:ndum),ny*ndum,MPI_REAL,right,tag,aay0(0,1:ny,1:ndum),ny*ndum,MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aay0(1:2,1:ny,1:ndum),2*ny*ndum,MPI_REAL,left,tag,aay0(nx+1:nx+2,1:ny,1:ndum),2*ny*ndum,MPI_REAL,right,tag,comm,stat,ierr)

       else
          call MPI_SENDRECV(aax0(nx,1:ny+1,1:ndum),(ny+1)*ndum,MPI_REAL,right,tag,aax0(0,1:ny+1,1:ndum),(ny+1)*ndum,MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aax0(1,1:ny+1,1:ndum),(ny+1)*ndum,MPI_REAL,left,tag,aax0(nx+1,1:ny+1,1:ndum),(ny+1)*ndum,MPI_REAL,right,tag,comm,stat,ierr)

          call MPI_SENDRECV(aay0(nx,1:ny,1:ndum),ny*ndum,MPI_REAL,right,tag,aay0(0,1:ny,1:ndum),ny*ndum,MPI_REAL,left,tag,comm,stat,ierr)
          call MPI_SENDRECV(aay0(1:2,1:ny,1:ndum),2*ny*ndum,MPI_REAL,left,tag,aay0(nx+1:nx+2,1:ny,1:ndum),2*ny*ndum,MPI_REAL,right,tag,comm,stat,ierr)

       endif

          !   do vertical transfer.

       if (nproc(2) .gt. 1) then

       call MPI_SENDRECV(aax0(0:nx+1,ny,1:ndum),(nx+2)*ndum,MPI_REAL,up,tag,aax0(0:nx+1,0,1:ndum),(nx+2)*ndum,MPI_REAL,down,tag,comm,stat,ierr)
       call MPI_SENDRECV(aax0(0:nx+1,1:2,1:ndum),(nx+2)*2*ndum,MPI_REAL,down,tag,aax0(0:nx+1,ny+1:ny+2,1:ndum),(nx+2)*2*ndum,MPI_REAL,up,tag,comm,stat,ierr)

       call MPI_SENDRECV(aay0(0:nx+2,ny,1:ndum),(nx+3)*ndum,MPI_REAL,up,tag,aay0(0:nx+2,0,1:ndum),(nx+3)*ndum,MPI_REAL,down,tag,comm,stat,ierr)
       call MPI_SENDRECV(aay0(0:nx+2,1,1:ndum),(nx+3)*ndum,MPI_REAL,down,tag,aay0(0:nx+2,ny+1,1:ndum),(nx+3)*ndum,MPI_REAL,up,tag,comm,stat,ierr)

       endif



       if (up .eq. MPI_PROC_NULL) then
            aax0(0:nx+1,ny+2,1:ndum) = 0.0
            aay0(0:nx+2,ny+1,1:ndum) = 0.0
         endif

       if (down .eq. MPI_PROC_NULL) then
            aax0(0:nx+1,0,1:ndum) = 0.0
            aay0(0:nx+2,0,1:ndum) = 0.0
       endif
       if (right .eq. MPI_PROC_NULL) then
            aax0(nx+1,0:ny+2,1:ndum) = 0.0
            aay0(nx+2,0:ny+1,1:ndum) = 0.0
         endif

       if (left .eq. MPI_PROC_NULL) then
            aax0(0,0:ny+2,1:ndum) = 0.0
            aay0(0,0:ny+1,1:ndum) = 0.0
       endif

endif

  end subroutine readdata

  subroutine writedata(nt,nr,tp)

    character*2 :: tp
    integer, intent(in) ::nt, nr
    integer :: i,j,k
    real, dimension(:,:,:), allocatable :: aax_global, aay_global, aaz_global
    integer, dimension(mpidir) :: dumcord


    if (tp .eq. '3d') then

       if (rank .eq. rankstart) then
          allocate(aax_global(nxglobal,nyglobal+1,nzglobal+1))
          allocate(aay_global(nxglobal+1,nyglobal,nzglobal+1))
          allocate(aaz_global(nxglobal+1,nyglobal+1,nzglobal))

          aax_global(1:nx,1:ny+1,1:nz+1) = aax(1:nx,1:ny+1,1:nz+1)
          aay_global(1:nx+1,1:ny,1:nz+1) = aay(1:nx+1,1:ny,1:nz+1)
          aaz_global(1:nx+1,1:ny+1,1:nz)   = aaz(1:nx+1,1:ny+1,1:nz)
       endif

       if (rank .ne. rankstart) then
            if ( coords(2) .ne. nproc(2)-1) then
                 if (coords(1) .ne. nproc(1)-1) then
                     CALL MPI_SEND(aax(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aay(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aaz(1:nx,1:ny,1:nz),nx*ny*nz,MPI_REAL,rankstart,tag,comm,ierr)
                 else
                     CALL MPI_SEND(aax(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aay(1:nx+1,1:ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aaz(1:nx+1,1:ny,1:nz),(nx+1)*ny*nz,MPI_REAL,rankstart,tag,comm,ierr)
                 endif
           else
                if (coords(1) .ne. nproc(1)-1) then
                     CALL MPI_SEND(aax(1:nx,1:ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aay(1:nx,1:ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aaz(1:nx,1:ny+1,1:nz),nx*(ny+1)*nz,MPI_REAL,rankstart,tag,comm,ierr)
                 else
                     CALL MPI_SEND(aax(1:nx,1:ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aay(1:nx+1,1:ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,rankstart,tag,comm,ierr)
                     CALL MPI_SEND(aaz(1:nx+1,1:ny+1,1:nz),(nx+1)*(ny+1)*nz,MPI_REAL,rankstart,tag,comm,ierr)
                 endif
           endif
     else

         do j=0,nproc(2)-1
             do i=0,nproc(1)-1
                dumcord(1) = i
                dumcord(2) = j
                CALL MPI_CART_RANK(comm,dumcord,nextrank,ierr)
                if (nextrank .ne. rankstart) then
                    if ( j .ne. nproc(2)-1) then
                        if (i .ne. nproc(1)-1) then
                             CALL MPI_RECV(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aay_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aaz_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz),nx*ny*nz,MPI_REAL,nextrank,tag,comm,stat,ierr)
                        else
                             CALL MPI_RECV(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aay_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aaz_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:nz),(nx+1)*ny*nz,MPI_REAL,nextrank,tag,comm,stat,ierr)
                        endif
                   else
                        if (i .ne. nproc(1)-1) then
                             CALL MPI_RECV(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aay_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny,1:nz+1),nx*ny*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aaz_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:nz),nx*(ny+1)*nz,MPI_REAL,nextrank,tag,comm,stat,ierr)
                        else
                             CALL MPI_RECV(aax_global((i*nx)+1:(i+1)*nx,(j*ny)+1:(j+1)*ny+1,1:nz+1),nx*(ny+1)*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aay_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny,1:nz+1),(nx+1)*ny*(nz+1),MPI_REAL,nextrank,tag,comm,stat,ierr)
                             CALL MPI_RECV(aaz_global((i*nx)+1:(i+1)*nx+1,(j*ny)+1:(j+1)*ny+1,1:nz),(nx+1)*(ny+1)*nz,MPI_REAL,nextrank,tag,comm,stat,ierr)
                        endif
                    endif
                endif
             enddo
         enddo
      endif

       if (rank .eq. rankstart) then

          write (filename,fmt='(a,"/",a,"_",i5.5,"p")')  &     ! *** Mac
              dir(1:length(dir,30)),root(1:length(root,10)),nt
          print *,'Saving model on file '//filename
         open(unit=1,file=filename,form='unformatted',status='unknown')
         write (1) opt
         write (1) (((aax_global(i,j,k),i=1,nxglobal  ),j=1,nyglobal+1),k=1,nzglobal+1)
         write (1) (((aay_global(i,j,k),i=1,nxglobal+1),j=1,nyglobal  ),k=1,nzglobal+1)
         write (1) (((aaz_global(i,j,k),i=1,nxglobal+1),j=1,nyglobal+1),k=1,nzglobal  )
         close(1)

         deallocate(aax_global,aay_global,aaz_global)
       endif
    endif

  end subroutine writedata

end module io
