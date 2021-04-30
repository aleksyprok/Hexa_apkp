module cal
   use var_global
  implicit none


  contains

!------------------------------------------------------------------
!  Compute magnetic field and electric current:
!------------------------------------------------------------------
  subroutine step!(mode1,numstep,base)

!  This routine takes one major step of the magneto-frictional relaxation.
!  It calculates the magnetic field B, current density C (=curl B),
!  and velocity V (= C x B / B**2), and then advances the vector
!  potential A according to an explicit scheme.

    use grid

    implicit none


!  Other variables:

    integer mode1,i,j,k,numstep,it
    real :: bbmax,eta,frc,local_max
    real :: localmin,globalmin
    character*1 :: base
    real :: t
    real :: b2local, b2global
    real, parameter :: c=3.e10 !cm/s
    real, dimension(nx+1,ny+1) :: exbase,eybase,bxbase,bybase,poynting
    real :: szlocal, szglobal
    real :: diss(nx+1,ny+1,nz+1), qlocal, qglobal
    real :: j2(nx+1,ny+1,nz+1),localj2,globalj2
    real :: localhyper,globalhyper
    real :: localj,globalj
    real :: alpha(nx+1,ny+1,nz+1)
    real :: linear(0:nx+2)


    frc=frc_coef*(1.E10)/(length_cm**2)*time_s !define frictional coefficient in dimensionless units

    t=0. !set local time variable to zero

    if (open .eq. 1) then !set up array for flux balancing
  	  do i=0,nx+2
  	    linear(i)=nx*coords(1)+(i-1) !effectively creates findgen(nx)
  	  enddo

  	  linear=linear*delx*imb !correction to Ay = imbalance/pixel*dx*x

  	  do k=1,nz+1
  	    do j=0,ny+1
  		  ay_corr(:,j,k)=linear(:) !Ay_corr = f(x) = imbalance/pixel*dx*x in sll volume
  	    enddo
  	  enddo

  	  daay=daay-ay_corr(:,:,1) !subtract correction from surface evolution ay (this is added on at the end to the whole volume)
    endif






do while (t .lt. 1.) ! loop until time has reached 1

    if (rank .eq. rankstart) print*,time,dt,basedt,dt/basedt

    t=t+dt
    time=time+dt


      aax(0:nx+1,0:ny+2,1)=aax(0:nx+1,0:ny+2,1)+daax(0:nx+1,0:ny+2)*dt
      aay(0:nx+2,0:ny+1,1)=aay(0:nx+2,0:ny+1,1)+daay(0:nx+2,0:ny+1)*dt


!  Compute div(A) at sides and corners


    div(1:nx+1,1:ny+1,2:nz)=(aax(1:nx+1,1:ny+1,2:nz)-aax(0:nx ,1:ny+1,2:nz  ))/delx  &
                       +(aay(1:nx+1,1:ny+1,2:nz)-aay(1:nx+1 ,0:ny ,2:nz  ))/dely  &
                       +(aaz(1:nx+1,1:ny+1,2:nz)-aaz(1:nx+1  ,1:ny+1  ,1:nz-1))/delz

    div(1:nx+1,1:ny+1,   1)=(aax(1:nx+1,1:ny+1,   1)-aax(0:nx,1:ny+1 ,     1))/delx  &
                         +(aay(1:nx+1,1:ny+1,   1)-aay(1:nx+1  ,0:ny ,     1))/dely

    if ( down .eq. MPI_PROC_NULL) then
      div(1:nx+1,1,1:nz+1) = 0.0
    endif
    if ( up .eq. MPI_PROC_NULL) then
      div(1:nx+1,ny+1,1:nz+1) = 0.0
    endif
    if ( right .eq. MPI_PROC_NULL) then
      div(nx+1,1:ny+1,1:nz+1) = 0.0
    endif
    if (left .eq. MPI_PROC_NULL) then
      div(1,1:ny+1,1:nz+1) = 0.0
    endif

    div(1:nx+1,1:ny+1,nz+1) = 0.0

!  Compute BBX, BBY, BBZ on cell faces:

    bbz(1:nx+2,1:ny+2,1:nz+1)=  &
       (aay(1:nx+2,0:ny+1,1:nz+1)-aay(0:nx+1,0:ny+1,1:nz+1))/delx  &
      -(aax(0:nx+1,1:ny+2,1:nz+1)-aax(0:nx+1,0:ny+1,1:nz+1))/dely

    if (down .eq. MPI_PROC_NULL) then
      bbz(1:nx+2,   1,1:nz+1)=bbz(1:nx+2,   2,1:nz+1)
    endif
    if (up .eq. MPI_PROC_NULL) then
      bbz(1:nx+2,ny+2,1:nz+1)=bbz(1:nx+2,ny+1,1:nz+1)
    endif
    if (left .eq. MPI_PROC_NULL) then
      bbz(1,1:ny+2,1:nz+1)=bbz(2,1:ny+2,1:nz+1)
    endif
    if (right .eq. MPI_PROC_NULL) then
      bbz(nx+2,1:ny+2,1:nz+1)=bbz(nx+1,1:ny+2,1:nz+1)
    endif



    bbx(1:nx+1,1:ny+2,2:nz+1)=  &
       (aaz(1:nx+1,1:ny+2,1:nz)-aaz(1:nx+1,0:ny+1,1:nz))/dely  &
      -(aay(1:nx+1,0:ny+1,2:nz+1)-aay(1:nx+1,0:ny+1,1:nz))/delz

    if (down .eq.MPI_PROC_NULL ) then
      bbx(1:nx+1,   1,2:nz+1)=bbx(1:nx+1,   2,2:nz+1)
    endif
    if (up .eq. MPI_PROC_NULL) then
      bbx(1:nx+1,ny+2,2:nz+1)=bbx(1:nx+1,ny+1,2:nz+1)
    endif

    bbx(1:nx+1,1:ny+2,   1)=bbx(1:nx+1,1:ny+2,   2)  &
         -delz/delx*(bbz(2:nx+2,1:ny+2,1)-bbz(1:nx+1,1:ny+2,1))

    if (open .eq. 1) then
        bbx(1:nx+1,1:ny+2,nz+2) = (bbz(2:nx+2,1:ny+2,nz+1) - bbz(1:nx+1,1:ny+2,nz+1))*delz/delx + bbx(1:nx+1,1:ny+2,nz+1)
    else
        bbx(1:nx+1,1:ny+2,nz+2)=bbx(1:nx+1,1:ny+2,nz+1)
    endif


    bby(1:nx+2,1:ny+1,2:nz+1)=  &
       (aax(0:nx+1,1:ny+1,2:nz+1)-aax(0:nx+1,1:ny+1,1:nz))/delz  &
      -(aaz(1:nx+2,1:ny+1,1:nz)-aaz(0:nx+1,1:ny+1,1:nz))/delx

    if (right .eq.MPI_PROC_NULL ) then
      bby(nx+2, 1:ny+1,2:nz+1)=bby(nx+1, 1:ny+1,2:nz+1)
    endif
    if (left .eq. MPI_PROC_NULL) then
      bby(1,1:ny+1,2:nz+1)=bby(2,1:ny+1,2:nz+1)
    endif

    bby(1:nx+2,1:ny+1,   1)=bby(1:nx+2,1:ny+1,   2)  &
         -delz/dely*(bbz(1:nx+2,2:ny+2,1)-bbz(1:nx+2,1:ny+1,1))
    if (open .eq. 1) then
       bby(1:nx+2,1:ny+1,nz+2) = (bbz(1:nx+2,2:ny+2,nz+1) - bbz(1:nx+2,1:ny+1,nz+1))*delz/dely + bby(1:nx+2,1:ny+1,nz+1)
    else
       bby(1:nx+2,1:ny+1,nz+2)=bby(1:nx+2,1:ny+1,nz+1)
    endif


!  Compute CCX, CCY, CCZ at edges:


    ccz(1:nx+1,1:ny+1,1:nz+2)=  &
         (bby(2:nx+2,1:ny+1,1:nz+2)-bby(1:nx+1,1:ny+1,1:nz+2))/delx  &
        -(bbx(1:nx+1,2:ny+2,1:nz+2)-bbx(1:nx+1,1:ny+1,1:nz+2))/dely

    ccx(1:nx+2,1:ny+1,1:nz+1)=  &
         (bbz(1:nx+2,2:ny+2,1:nz+1)-bbz(1:nx+2,1:ny+1,1:nz+1))/dely  &
        -(bby(1:nx+2,1:ny+1,2:nz+2)-bby(1:nx+2,1:ny+1,1:nz+1))/delz

    ccy(1:nx+1,1:ny+2,1:nz+1)=  &
         (bbx(1:nx+1,1:ny+2,2:nz+2)-bbx(1:nx+1,1:ny+2,1:nz+1))/delz  &
        -(bbz(2:nx+2,1:ny+2,1:nz+1)-bbz(1:nx+1,1:ny+2,1:nz+1))/delx

!  Magnetic field at cell corners:

    bx(1:nx+1,1:ny+1,1:nz+1)=  &
         0.25*(bbx(1:nx+1,1:ny+1,1:nz+1)+bbx(1:nx+1,2:ny+2,1:nz+1)  &
              +bbx(1:nx+1,1:ny+1,2:nz+2)+bbx(1:nx+1,2:ny+2,2:nz+2))

    by(1:nx+1,1:ny+1,1:nz+1)=  &
         0.25*(bby(1:nx+1,1:ny+1,1:nz+1)+bby(2:nx+2,1:ny+1,1:nz+1)  &
              +bby(1:nx+1,1:ny+1,2:nz+2)+bby(2:nx+2,1:ny+1,2:nz+2))

    bz(1:nx+1,1:ny+1,1:nz+1)=  &
         0.25*(bbz(1:nx+1,1:ny+1,1:nz+1)+bbz(2:nx+2,1:ny+1,1:nz+1)  &
              +bbz(1:nx+1,2:ny+2,1:nz+1)+bbz(2:nx+2,2:ny+2,1:nz+1))

    bb=bx*bx+by*by+bz*bz



!  Field strength with minimum value for magneto-friction:

    do k=1,nz+1
       local_max=0.0
       do j=1,ny+1 !changed order for more efficient memory access
          do i=1,nx+1 !changed order for more efficient memory access
             local_max=max(bb(i,j,k),local_max)
          enddo
       enddo
       CALL MPI_ALLREDUCE(local_max,bbmax,1,MPI_REAL,MPI_MAX,comm,ierr)
       bbmax=1.e-4*bbmax
       do j=1,ny+1
          do i=1,nx+1
             bbm(i,j,k)=max(bb(i,j,k),bbmax)
          enddo
       enddo
    enddo

!  Current density at cell corners:
    cx(1:nx+1,1:ny+1,1:nz+1)=  &
         0.5*(ccx(1:nx+1,1:ny+1,1:nz+1)+ccx(2:nx+2,1:ny+1,1:nz+1))

    cy(1:nx+1,1:ny+1,1:nz+1)=  &
         0.5*(ccy(1:nx+1,1:ny+1,1:nz+1)+ccy(1:nx+1,2:ny+2,1:nz+1))

    cz(1:nx+1,1:ny+1,1:nz+1)=  &
         0.5*(ccz(1:nx+1,1:ny+1,1:nz+1)+ccz(1:nx+1,1:ny+1,2:nz+2))

    ch=(bx*cx+by*cy+bz*cz)/bbm

    j2=cx*cx+cy*cy+cz*cz

    alpha=ch

!  Apply isotropic diffusion of div(A):

    if(etad.gt.0.0) then
       eta=etad*dt

       if ( coords(2) .ne. nproc(2)-1) then
          aax(1:nx  ,1:ny  ,1:nz)=aax(1:nx  ,1:ny  ,1:nz)  &
              +eta*(div(2:nx+1,1:ny  ,1:nz )-div(1:nx  ,1:ny  ,1:nz))/delx

          if  ( coords(1) .ne. nproc(1)-1) then
            aay(1:nx,1:ny  ,1:nz)=aay(1:nx,1:ny  ,1:nz)  &
                 +eta*(div(1:nx,2:ny+1,1:nz  )-div(1:nx,1:ny  ,1:nz))/dely
            aaz(1:nx,1:ny,1:nz)=aaz(1:nx,1:ny,1:nz)  &
                 +eta*(div(1:nx,1:ny,2:nz+1)-div(1:nx,1:ny,1:nz))/delz
          else
            aay(1:nx+1,1:ny  ,1:nz)=aay(1:nx+1,1:ny  ,1:nz)  &
                 +eta*(div(1:nx+1,2:ny+1,1:nz  )-div(1:nx+1,1:ny  ,1:nz))/dely
            aaz(1:nx+1,1:ny,1:nz)=aaz(1:nx+1,1:ny,1:nz)  &
                 +eta*(div(1:nx+1,1:ny,2:nz+1)-div(1:nx+1,1:ny,1:nz))/delz
          endif
       else
           aax(1:nx  ,1:ny+1  ,1:nz)=aax(1:nx  ,1:ny+1  ,1:nz)  &
              +eta*(div(2:nx+1,1:ny+1 ,1:nz )-div(1:nx  ,1:ny+1  ,1:nz))/delx

          if  ( coords(1) .ne. nproc(1)-1) then
            aay(1:nx,1:ny  ,1:nz)=aay(1:nx,1:ny  ,1:nz)  &
                 +eta*(div(1:nx,2:ny+1,1:nz)-div(1:nx,1:ny  ,1:nz))/dely
            aaz(1:nx,1:ny+1,1:nz)=aaz(1:nx,1:ny+1,1:nz)  &
                 +eta*(div(1:nx,1:ny+1,2:nz+1)-div(1:nx,1:ny+1,1:nz))/delz
          else
            aay(1:nx+1,1:ny  ,1:nz)=aay(1:nx+1,1:ny  ,1:nz)  &
                 +eta*(div(1:nx+1,2:ny+1,1:nz  )-div(1:nx+1,1:ny  ,1:nz))/dely
            aaz(1:nx+1,1:ny+1,1:nz)=aaz(1:nx+1,1:ny+1,1:nz)  &
                 +eta*(div(1:nx+1,1:ny+1,2:nz+1)-div(1:nx+1,1:ny+1,1:nz))/delz
          endif
        endif


    endif

!  Apply isotropic diffusion (not at lower/upper boundaries):

    if(etai.gt.0.0) then
       eta=etai*dt

      if ( coords(2) .ne. nproc(2)-1) then
         aax(1:nx  ,1:ny,2:nz)=aax(1:nx  ,1:ny,2:nz)  &
                          -eta*ccx(2:nx+1,1:ny,2:nz)

         if  ( coords(1) .ne. nproc(1)-1) then
           aay(1:nx,1:ny  ,2:nz)=aay(1:nx,1:ny  ,2:nz)  &
                          -eta*ccy(1:nx,2:ny+1,2:nz)
           aaz(1:nx,1:ny,1:nz)=aaz(1:nx,1:ny,1:nz)  &
                          -eta*ccz(1:nx,1:ny,2:nz+1)
         else
           aay(1:nx+1,1:ny  ,2:nz)=aay(1:nx+1,1:ny  ,2:nz)  &
                          -eta*ccy(1:nx+1,2:ny+1,2:nz)
           aaz(1:nx+1,1:ny,1:nz)=aaz(1:nx+1,1:ny,1:nz)  &
                          -eta*ccz(1:nx+1,1:ny,2:nz+1)
         endif

      else

         aax(1:nx  ,1:ny+1,2:nz)=aax(1:nx  ,1:ny+1,2:nz)  &
                          -eta*ccx(2:nx+1,1:ny+1,2:nz)

         if  ( coords(1) .ne. nproc(1)-1) then
           aay(1:nx,1:ny  ,2:nz)=aay(1:nx,1:ny  ,2:nz)  &
                          -eta*ccy(1:nx,2:ny+1,2:nz)
           aaz(1:nx,1:ny+1,1:nz)=aaz(1:nx,1:ny+1,1:nz)  &
                          -eta*ccz(1:nx,1:ny+1,2:nz+1)
         else
           aay(1:nx+1,1:ny  ,2:nz)=aay(1:nx+1,1:ny  ,2:nz)  &
                          -eta*ccy(1:nx+1,2:ny+1,2:nz)
           aaz(1:nx+1,1:ny+1,1:nz)=aaz(1:nx+1,1:ny+1,1:nz)  &
                          -eta*ccz(1:nx+1,1:ny+1,2:nz+1)
         endif

       endif

    endif

!  Magneto-frictional velocity:

    vx=frc*(cy*bz-cz*by)/bbm
    vy=frc*(cz*bx-cx*bz)/bbm
    vz=frc*(cx*by-cy*bx)/bbm

!  Rate of change of vector potential:

    cx=vy*bz-vz*by
    cy=vz*bx-vx*bz
    cz=vx*by-vy*bx

!  Helicity fluxes, and effect of hyperdiffusion:

    if(eta4.gt.0.0) then
    !B^2 grad alpha
       ccx(2:nx+1,1:ny+1,1:nz+1)=  &
            0.5*(bb(2:nx+1,1:ny+1,1:nz+1)+bb(1:nx,1:ny+1,1:nz+1))  &
               *(ch(2:nx+1,1:ny+1,1:nz+1)-ch(1:nx,1:ny+1,1:nz+1))/delx

    call MPI_SENDRECV(ccx(nx+1,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,right,tag,ccx(1,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,left,tag,comm,stat,ierr)
    call MPI_SENDRECV(ccx(2,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,left,tag,ccx(nx+2,1:ny+1,1:nz+1),(ny+1)*(nz+1),MPI_REAL,right,tag,comm,stat,ierr)

       if (right .eq. MPI_PROC_NULL ) then
         ccx(nx+2,1:ny+1 ,1:nz+1)=ccx(nx+1, 1:ny+1,1:nz+1)
       endif
       if (left .eq. MPI_PROC_NULL ) then
        ccx(1,1:ny+1,1:nz+1)=ccy(2,1:ny+1,1:nz+1)
       endif

       ccy(1:nx+1,2:ny+1,1:nz+1)=  &
            0.5*(bb(1:nx+1,2:ny+1,1:nz+1)+bb(1:nx+1,1:ny,1:nz+1))  &
               *(ch(1:nx+1,2:ny+1,1:nz+1)-ch(1:nx+1,1:ny,1:nz+1))/dely

    call MPI_SENDRECV(ccy(1:nx+1,ny+1,1:nz+1),(nx+1)*(nz+1),MPI_REAL,up,tag,ccy(1:nx+1,1,1:nz+1),(nx+1)*(nz+1),MPI_REAL,down,tag,comm,stat,ierr)
    call MPI_SENDRECV(ccy(1:nx+1,2,1:nz+1),(nx+1)*(nz+1),MPI_REAL,down,tag,ccy(1:nx+1,ny+2,1:nz+1),(nx+1)*(nz+1),MPI_REAL,up,tag,comm,stat,ierr)

       if (down .eq. MPI_PROC_NULL ) then
         ccy(1:nx+1,   1,1:nz+1)=ccy(1:nx+1,   2,1:nz+1)
       endif
       if ( up .eq. MPI_PROC_NULL ) then
        ccy(1:nx+1,ny+2,1:nz+1)=ccy(1:nx+1,ny+1,1:nz+1)
       endif

       ccz(1:nx+1,1:ny+1,2:nz+1)=  &
            0.5*(bb(1:nx+1,1:ny+1,2:nz+1)+bb(1:nx+1,1:ny+1,1:nz))  &
               *(ch(1:nx+1,1:ny+1,2:nz+1)-ch(1:nx+1,1:ny+1,1:nz))/delz

       ccz(1:nx+1,1:ny+1,   1)=ccz(1:nx+1,1:ny+1,   2)
       ccz(1:nx+1,1:ny+1,nz+2)=ccz(1:nx+1,1:ny+1,nz+1)

       ch=eta4*(  (ccx(2:nx+2,1:ny+1,1:nz+1)-ccx(1:nx+1,1:ny+1,1:nz+1))/delx  &
               +  (ccy(1:nx+1,2:ny+2,1:nz+1)-ccy(1:nx+1,1:ny+1,1:nz+1))/dely  &
               +  (ccz(1:nx+1,1:ny+1,2:nz+2)-ccz(1:nx+1,1:ny+1,1:nz+1))/delz) &
               /bbm

       if (down .eq. MPI_PROC_NULL ) then
         ch(1:nx+1,1,1:nz+1) = 0.0
       endif
       if (up .eq. MPI_PROC_NULL ) then
         ch(1:nx+1,ny+1,1:nz+1) = 0.0
       endif
       if (right .eq. MPI_PROC_NULL ) then
         ch(nx+1,1:ny+1,1:nz+1) = 0.0
       endif
       if (left .eq. MPI_PROC_NULL ) then
         ch(1,1:ny+1,1:nz+1) = 0.0
       endif

         ch(1:nx+1,1:ny+1,nz+1) = 0.0

       cx=cx+ch*bx
       cy=cy+ch*by
       cz=cz+ch*bz

    endif

!  Apply advection and hyper-diffusion (unless mode=1):



    if(mode1.ne.1) then

       if ( coords(2) .ne. nproc(2)-1) then

           aax(1:nx  ,1:ny,2:nz+1)=aax(1:nx  ,1:ny,2:nz+1)  &
            +0.5*dt*(cx(2:nx+1,1:ny,2:nz+1)+cx(1:nx,1:ny,2:nz+1))

           if ( coords(1) .ne. nproc(1)-1) then
             aay(1:nx,1:ny  ,2:nz+1)=aay(1:nx,1:ny  ,2:nz+1)  &
                 +0.5*dt*(cy(1:nx,2:ny+1,2:nz+1)+cy(1:nx,1:ny,2:nz+1))
             aaz(1:nx,1:ny,1:nz)=aaz(1:nx,1:ny,1:nz)  &
                 +0.5*dt*(cz(1:nx,1:ny,2:nz+1)+cz(1:nx,1:ny,1:nz))
           else
             aay(1:nx+1,1:ny  ,2:nz+1)=aay(1:nx+1,1:ny  ,2:nz+1)  &
                 +0.5*dt*(cy(1:nx+1,2:ny+1,2:nz+1)+cy(1:nx+1,1:ny,2:nz+1))
             aaz(1:nx+1,1:ny,1:nz)=aaz(1:nx+1,1:ny,1:nz)  &
                 +0.5*dt*(cz(1:nx+1,1:ny,2:nz+1)+cz(1:nx+1,1:ny,1:nz))
           endif
       else
           aax(1:nx  ,1:ny+1,2:nz+1)=aax(1:nx  ,1:ny+1,2:nz+1)  &
            +0.5*dt*(cx(2:nx+1,1:ny+1,2:nz+1)+cx(1:nx,1:ny+1,2:nz+1))

           if ( coords(1) .ne. nproc(1)-1) then
             aay(1:nx,1:ny  ,2:nz+1)=aay(1:nx,1:ny  ,2:nz+1)  &
                 +0.5*dt*(cy(1:nx,2:ny+1,2:nz+1)+cy(1:nx,1:ny,2:nz+1))
             aaz(1:nx,1:ny+1,1:nz)=aaz(1:nx,1:ny+1,1:nz)  &
                 +0.5*dt*(cz(1:nx,1:ny+1,2:nz+1)+cz(1:nx,1:ny+1,1:nz))
           else
             aay(1:nx+1,1:ny  ,2:nz+1)=aay(1:nx+1,1:ny  ,2:nz+1)  &
                 +0.5*dt*(cy(1:nx+1,2:ny+1,2:nz+1)+cy(1:nx+1,1:ny,2:nz+1))
             aaz(1:nx+1,1:ny+1,1:nz)=aaz(1:nx+1,1:ny+1,1:nz)  &
                 +0.5*dt*(cz(1:nx+1,1:ny+1,2:nz+1)+cz(1:nx+1,1:ny+1,1:nz))
           endif
        endif
    endif

   if (open .eq. 1) then
     aay=aay+ay_corr*dt !correct Ay throughout the whole volume to ensure flux balance
   endif


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
            aax(0,1:ny+2,1:nz+1) = 0.0
            aay(0,1:ny+1,1:nz+1) = 0.0
            aaz(0,1:ny+2,1:nz) = 0.0
       endif

       !correction to Ay ghost cells if periodic and open BCs are used.
  if (periodic .eq. 1 .and. open .eq. 1) then
	if (coords(1) .eq. 0) then !left most processors
	  aay(0,:,:)=aay(0,:,:)-nxglobal*delx*(imb*t+totimb)
	endif

	if (coords(1) .eq. nproc(1)-1) then !right most processors
	  aay(nx+1,:,:)=aay(nx+1,:,:) +nxglobal*delx*(imb*t+totimb)
	  aay(nx+2,:,:)=aay(nx+2,:,:) +nxglobal*delx*(imb*t+totimb)
	endif
  endif

!########################################################
!Calculate and write diagnostic data to file
!########################################################


      ! calculate magnetic energy

      b2local=sum(bb(1:nx,1:ny,1:nz+1))/8./pi*delx**3

      CALL MPI_ALLREDUCE(b2local,b2global,1,MPI_Real,MPI_SUM,MPI_COMM_WORLD,ierr)




      !poynting flux c/4pi * (ExBy-EyBx) where E=-dA/dt (or -v x B)

      bxbase=bx(:,:,2)
      bybase=by(:,:,2)

      exbase=-0.5*(daax(0:nx,1:ny+1)+daax(1:nx+1,1:ny+1))
      eybase=-0.5*(daay(1:nx+1,0:ny)+daay(1:nx+1,1:ny+1))

      exbase=-1.*(vy(:,:,2)*bz(:,:,2)  -vz(:,:,2)*by(:,:,2))
      eybase=-1.*(vz(:,:,2)*bx(:,:,2)  -vx(:,:,2)*bz(:,:,2))

      poynting = 0.

      poynting = c/4./pi * (exbase*bybase-eybase*bxbase)*delx**2
      szlocal=sum(poynting(1:nx,1:ny))

      CALL MPI_ALLREDUCE(szlocal,szglobal,1,MPI_Real,MPI_SUM,MPI_COMM_WORLD,ierr)




      !dissipation Q=B^2/4pi * nu * v^2    !nu=1/frc

      diss=vx*vx + vy*vy + vz*vz
      diss=diss/frc/4./pi
      diss=diss*bb

      qlocal=sum(diss(1:nx,1:ny,1:nz+1))*delx**3
      CALL MPI_ALLREDUCE(qlocal,qglobal,1,MPI_Real,MPI_SUM,MPI_COMM_WORLD,ierr)




      !Ohmic heating
      localj2=sum(j2(1:nx,1:ny,1:nz+1))*etai*delx**3/4./pi
      CALL MPI_ALLREDUCE(localj2,globalj2,1,MPI_Real,MPI_SUM,MPI_COMM_WORLD,ierr)



      !hyperdiffusion Q = B^2/4/pi eta4 |grad(alpha)|^2
      if (eta4 .gt. 0) then
	!grad alpha (x-ccx, y-ccy,z-ccz)
	ccx(2:nx+1,:,:) = (alpha(2:nx+1,:,:)-alpha(1:nx,:,:))/delx
	ccx(1,:,:) = ccx(2,:,:)
	ccx(nx+2,:,:) = ccx(nx+1,:,:)

	ccy(:,2:ny+1,:) = (alpha(:,2:ny+1,:)-alpha(:,1:ny,:))/dely
	ccy(:,1,:) = ccy(:,2,:)
	ccx(:,ny+2,:) = ccx(:,ny+1,:)

	ccz(:,:,2:nz+1) = (alpha(:,:,2:nz+1)-alpha(:,:,1:nz))/delz
	ccz(:,:,1) = ccz(:,:,2)
	ccz(:,:,nz+2) = ccz(:,:,nz+1)

	!|grad alpha|^2
	alpha(:,:,:) = (ccx(1:nx+1,:,:)+ccx(2:nx+2,:,:))**2/4. &
		    + (ccy(:,1:ny+1,:)+ccy(:,2:ny+1,:))**2/4. &
		    + (ccz(:,:,1:nz+1)+ccz(:,:,2:nz+2))**2/4.

	alpha = alpha*bb*eta4


	localhyper=sum(abs(alpha(1:nx,1:ny,1:nz+1)))*delx**3/4./pi
      else
        localhyper=0.
      endif
      CALL MPI_ALLREDUCE(localhyper,globalhyper,1,MPI_Real,MPI_SUM,MPI_COMM_WORLD,ierr)



      !current
      localj=sqrt(sum(j2(1:nx,1:ny,1:nz+1)))*delx**3
      CALL MPI_ALLREDUCE(localj,globalj,1,MPI_Real,MPI_SUM,MPI_COMM_WORLD,ierr)


      if (rank .eq. rankstart) then
	    write(50,*) time, & !time
	    dt*time_s,&  !timestep(s)
	    b2global*(length_cm**3),& !magnetic energy (erg)
	    szglobal*(length_cm**2)/time_s,& !poynting flux (erg/s)
	    qglobal/time_s*length_cm**3,& !frictional dissipation (erg/s)
	    globalj2*length_cm**3/time_s,& !Ohmic heating
	    globalhyper*length_cm**3/time_s,& !hyperdiffusion
	    globalj*length_cm**2 !currents
      endif







      !####################################################
      !                 Determine timestep
      !####################################################

      !determine minimum cell crossing time for advection and diffusion terms

      localmin=minval( (/ delx/maxval(vx),&
			  dely/maxval(vy),&
			  delz/maxval(vz),&
			  delx**2/etai,&
			  delx**2/etad,&
			  delx**4/eta4 /) )

      CALL MPI_ALLREDUCE(localMin,globalMin,1,MPI_Real,MPI_MIN,MPI_COMM_WORLD,ierr)

      !set timestep

	dt=globalmin*0.2

	if (dt .lt. basedt*1.e-3) dt=1.e-3*basedt !if timestep is too small, set to 0.001*basedt
	if (t .gt. 1.) exit
	if (1.-t .lt. dt) dt=1.-t+1.e-3*basedt !if timestep is large enough to increase t above 1 then set dt to make t=1




enddo




end subroutine step

end module cal
