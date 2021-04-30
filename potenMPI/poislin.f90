      subroutine poislin
! *
! ***************************************************************
! *                                                             *
! *       MULTI-GRID poisson solver - Linear & 2D               *
! *                                                             *
! *  Solves the flux function equation:                         *
! *                                                             *
! *      d2(psi)/dx^2 + d2(psi)/dy^2 = -Bz(x,y)                 *
! *                                                             *
! *  by Guass-Seidel & Linear Multigrid with V cycles.          *
! *  Then calculates Ax & Ay  (A =curl(psi Z) ) and dumps them  *
! *  to file to be used as the base boundary conditions in the  *
! *  3D linear force-free field multigrid solver "f3cray.f".    *
! *                                                             *
! *       Ax = d(psi)/dy   Ay = -d(psi)/dx                      *
! *                                                             *
! *    -- Aaron Longbottom 25th July 1996                       *
! *                                                             *
! *    Updated                                                  *
! *    -- Duncan H Mackay   9th May 97                          *
! ***************************************************************
! *
! *******This line is 72 columns wide*************************************
! *
  use sysdat_mod
  use const_mod
  use xb_mod
  use xn_mod
  use times_mod
  use output_mod
  use rhs1_mod
  use ina_mod
  use wu_mod
  use res_mod
  use acc_mod
 
 implicit double precision(a-h,o-z)

! *
! **  Define constants used
! *
       ina=1
! *
       pi = 2.d0*dacos(0.d0)
       nx=na
       dx=1.0d0/dfloat(nx)
! *
! **  Input parameters for this run 
! *
       ncycle = 30     ! Maximum number of Multigrid cycles.
       acc    = 1.d-10 ! Required accuracy.
       
       
! ** Set up rhs of equation on top grid
       do 20 i=-nx,nx
       do 20 j=-nx,nx
          call mfinder(nx,i,j,m)
          ff(m)=0.d0
          xx = dfloat(i)*dx
          yy = dfloat(j)*dx
          xb(m) = (dcos(0.5d0*pi*xx))*(dcos(0.5d0*pi*yy))
20     continue
!*
       call bzfinder() !Reads magnetogram into ff

! **  Replace xn with xb
! *
       call replace(xn,xb,nm)

! **  Output current variables
! *
       iocall = 0
       tlast = time
       wu=0.d0


! *
       do 60 icoo = 0 ,ncycle ! outer loop
! *
          xiout  = dfloat(icoo)/dfloat(ncycle)
! *
! ** MAIN LOOP for set up begins here
! *
! ***********************************************************
! ** MULTIGRID CYCLE
! *
! ** Multigrid cycle
! *
           call driver()
! *
! ** Normal step
! *
! c          call replace(xb,xn,nm)
! c          call smoother(ff(1),xb(1),xn(1),res(1),nm,nx,dx)
! c          wu = wu + 1.d0
! *
! ***********************************************************
! *
          resmax=-1.d0
          do 2 m=1,nm
             pp=dabs(res(m))
             if (pp.gt.resmax) resmax=pp
2         continue
! *
          per=0.d0
          do 3 i=1,nm
             per=per+res(i)**2
3         continue
          per=(per/dfloat(nm))**0.5d0
! *
          time = time + dt
! *
! ***********************************************************
! *
! ** Write current values to file.
! *
          tlast = time

          call flush(nout)
! *
! ** Accuracy check
! *
          if (per.lt.acc) goto 70
! *
998       format(7x,'WU',10x,'Max Res',7x,'|Res|')
999       format(3(2x,d12.6))
! *
60     continue   !  Update compleated
! *
! ** Final calculation of outputs and dump.
! *
70     continue

       call adump()


       return
       end ! poislin.
       
       
       
! *
! *
! *******This line is 72 columns wide*************************************
! *
       subroutine adump()
! *
! ***********************************
! *  Calculates Ax and Ay and dumps *
! * them to binary file ready to be *
! * read by f3cray.f.               *
! ***********************************
! *
       use sysdat_mod 
       use avalues_mod
       use const_mod
       use xn_mod
       implicit double precision(a-h,o-z)


       tdx = 2.d0*dx
! *
      ! do 10 i = -nx,nx
      ! do 10 j = -nx,nx

! *
       do 10 i=1,2*nx+1
       do 10 j=1,2*nx+1
       
! ** Ay = -d(psi)/dx
! *
          if (i-nx-1.eq.-nx.or.i-nx-1.eq.nx) then
             ay(i,j) = 0.d0
          else
             call mfinder(nx,i+1-nx-1,j-nx-1,mp)
             mm = mp-2
             ay(i,j) = -(xn(mp)-xn(mm))/tdx
          endif
! *
! ** Ax = d(psi)/dy
! *
          if (j-nx-1.eq.-nx.or.j-nx-1.eq.nx) then
             ax(i,j) = 0.d0
          else
            call mfinder(nx,i-nx-1,j+1-nx-1,mp)
             mm = mp-4*nx-2
             ax(i,j) = (xn(mp)-xn(mm))/tdx
          endif
10     continue

! *
! ** Dump Ax & Ay to file.
! *
      ! open(unit=2,file='axpot.dat')
      ! write (2,99) ((ax(i,j),i=-nx,nx),j=-nx,nx)
      ! close(2)
! *
      ! open(unit=2,file='aypot.dat')
      ! write (2,99) ((ay(i,j),i=-nx,nx),j=-nx,nx)
      ! close(2)
! *
 99    format(1pe21.14)
! *
       return
       end ! adump.
! *
! *******This line is 72 columns wide**************************************
       subroutine replace(yn,yp,lyn)
! *
! *********************************
! **  Replaces array yn with yp  **
! *********************************
! *
       implicit double precision(a-h,o-z)
! *
       dimension yp(lyn),yn(lyn)
! *
       do 10 m = 1,lyn
          yn(m) = yp(m)
10     continue
!*
       return
       end ! replace


      SUBROUTINE ZEROARR(N,A)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION A(N)
      INTEGER I
      DO I=1,N
        A(I)=0
      ENDDO
      END

      SUBROUTINE BASE(NX,NY,NZ,RHSZ,AX,AY)
      IMPLICIT NONE
      INTEGER NX,NY,NZ
      DOUBLE PRECISION AX(NX,NY,NZ),AY(NX,NY,NZ),RHSZ(NX,NY,NZ)
      INTEGER I,J,K
      K=1
      DO J=2,NY-1
        DO I=2,NX-1
          RHSZ(I,J,K)=-(AX(I+1,J,K)-AX(I-1,J,K))-(AY(I,J+1,K)-AY(I,J-1,K))
        ENDDO
      ENDDO
      END
      
      
      
      
      
      
             subroutine driver()
!*
!*****************************************************
!**                                                 **
!**  Driver for multigrid sweep using V-cycle       **
!**  For solving LINEAR problem.                    **
!**                                                 **
!*****************************************************
       use sysdat_mod
       use ina_mod
       implicit double precision(a-h,o-z)

       integer          nk(nmax) 

!*
!** Set up initial values
!*
       igamma   = 1     ! V cycle
       nk(nmax) = 1     ! number of cycles
       kk       = nmax  ! top level
       k        = kk    ! current level
       ina      = 1
!*
10     if (nk(kk).gt.0) then
!*
!** Start F cycle
!*
          if (nk(k).eq.0.or.k.eq.1) then
!*
             if (k.eq.1) then
                call Sroutine(k,nk)
             endif
!*
             if (k.eq.kk) then
             else
!*
!** Goto finer grid
!*
                k = k+1
                call Broutine(k,nk)
             endif
             nk(k)=nk(k)-1
          else
!*
!** Goto coarser grid
!*
             call Aroutine(k,nk)
!*
             k = k-1
             nk(k) = igamma
          endif
!*
          goto 10
       endif
!*
!c      print*,'*************************************'
!*
       return
       end ! driver.
!*                
!**************** This is 72 columns wide ******************************* 
!*
       subroutine Aroutine(k,nk)
!*
        use sysdat_mod
        use const_mod
        use rhs1_mod
        use s_mod
        use rhs_mod
        use xcb_mod
        use xb_mod
        use xn_mod
        use times_mod
        use ina_mod
        use wu_mod
        use res_mod
        
       implicit double precision(a-h,o-z)

       double precision dthold(nmax)
       integer          nk(nmax)
!
!*
!** Smooth. Move down grid.
!*
!*
!** If on top grid use top grid smoother & calculate residue
!*
       if (k.eq.nmax) then  ! On top grid
!*
!** Zero xct for this cycle
!*
          do 1 m=1,nm/2
             xct(m)=0.d0
1         continue
!*
          ina=1
          call key(k-1,nmax,nst1,ncx1,lb1)
          call replace(xb(1),xn(1),nm)
          dcx=dx
          call smoother(ff(1),xb(1),xn(1),res(1),nm,nx,dx)
          wu=wu+2.d0**(2.d0*dfloat(k-nmax))
!*
!** Else use lower grid smoother & calculate residue
!*
       else                 ! On lower grid
          call key(k,nmax,nst,ncx,lb)
          call key(k-1,nmax,nst1,ncx1,lb1)
!*
          if (ina.eq.2) then
             call replace(xcb(nst),xcn(nst),lb)
          else
             call replace(xcb(nst),xct(nst),lb)
          endif
          ina=1
          dcx=1.0d0/dfloat(ncx)
!*
          call smoother(ffc(nst),xcb(nst),xcn(nst),resc(nst),&
                          lb,ncx,dcx)
          wu=wu+2.d0**(2.d0*dfloat(k-nmax)) 
!*
       endif
!*
!** Restrict solution just calculated & residue
!*
       resmax =-1.d0
       resmax1=-1.d0
       if (k.eq.nmax) then  ! On top grid
          call restrict(res(1),nm,resc(nst1),lb1,ncx1)
       else                 ! On lower grid
          call restrict(resc(nst),lb,resc(nst1),lb1,ncx1)
       endif
!*
!** Calculate f (ffc) on courser grid
!*
       dcx=1.0d0/dfloat(ncx1)
!*
       call fcal(1.d0,ffc(nst1),resc(nst1),xct(nst1),lb1,&
                ncx1,dcx)
!*
!*
!*
       return
       end! Aroutine.
!*
!**************** This is 72 columns wide ******************************* 
!*
       subroutine Broutine(k,nk) 
       
        use sysdat_mod
        use const_mod
        use rhs1_mod
        use rhs_mod
        use xcb_mod
        use xb_mod
        use xn_mod
        use times_mod
        use ina_mod
        use wu_mod
        use res_mod
!*
       implicit double precision(a-h,o-z)

       double precision s(nmax)
       integer          nk(nmax) 

       ina=2
!*
!** Move up grid. Smooth.
!* 
!*
!** Difference solution on this grid [ie. find u(n+1) - u(n)]
!*
       call key(k-1,nmax,nst,ncx,lb)
       do 10 i = nst,nst+lb-1
          xcb(i) = xcn(i)                   ! for linear problems
 10    continue
!*
!** Prolongate difference up to finer grid and add it to current u
!*
       if (k.eq.nmax) then ! On second top grid
          call prolongate(xb(1),nm,xcb(nst),lb,ncx)
!*
          do 20 j = 1,nm
             xb(j) = xb(j) + xn(j)
 20       continue
       else                  ! On lower grid
          call key(k,nmax,nst1,ncx1,lb1)
          call prolongate(xcb(nst1),lb1,xcb(nst),lb,ncx)
!*
          do 30 j = nst1,lb1+nst1-1
             xcb(j) = xcb(j) + xcn(j) 
 30       continue 
       endif
!* 
!** If on top grid calculate new rhs of equation
!** and use top grid smoother
!*
       if (k.eq.nmax) then
          call smoother(ff(1),xb(1),xn(1),res(1),nm,nx,dx)
          wu=wu+2.d0**(2.d0*dfloat(k-nmax)) 
!*
!** If on lower grid use lower grid smoother
!*
       else
          dcx=1.0d0/dfloat(ncx1)
          call smoother(ffc(nst1),xcb(nst1),xcn(nst1),&
                       resc(nst1),lb1,ncx1,dcx)
          wu=wu+2.d0**(2.d0*dfloat(k-nmax)) 
       endif
!*
!*
       return 
       end! Broutine. 
!*
!**************** This is 72 columns wide *******************************
!*
       subroutine Sroutine(k,nk) 
       
        use sysdat_mod
        use const_mod
        use acc_mod
        use rhs_mod
        use xcb_mod
        use xb_mod
        use xn_mod
        use times_mod
        use output_mod
        use ina_mod
        use wu_mod
        use res_mod
!*
       implicit double precision(a-h,o-z)

       integer          nk(nmax) 
!*
!** Smooth on grid till gain accuracy required.
!* 
       ina=3
!c      TOLER=1.d-10 !TOLERENCE WANTED ON COURSEST GRID
       ITMAX=500    !NEED TO CHANGE THIS FOR A REAL RUN
!*
!c      print*,'_____________________'
!*
!** Loop using lower grid smoother till required accuracy obtained.
!*
       call key(k,nmax,nst,ncx,lb)
       dcx=1.0d0/dfloat(ncx)
       dthold=dt
       dt=0.25d0*dcx*dcx
!*
       call replace(xcn(nst),xct(nst),lb)
!*
       do 10 n=1,itmax
!*
          call replace(xcb(nst),xcn(nst),lb)
          call smoother(ffc(nst),xcb(nst),xcn(nst),resc(nst),&
                          lb,ncx,dcx)
          wu=wu+2.d0**(2.d0*dfloat(k-nmax)) 
!*
!** Error (tolerence) check on current solution
!*
!c         resmax =-1.d0
!c         do 15 m=nst,nst+lb-1
!c            pp=dabs(resc(m))
!c            if (pp.gt.resmax ) resmax =pp
!c15       continue
!*
          per=0.d0
          do 5 i=nst,nst+lb-1
             per=per+resc(i)**2
 5        continue 
          per=(per/dfloat(lb))**0.5d0
!*
!************************************************************ 
!* 
!c         print*,per,resmax
!c         print*,'S ',k,nk(k),wu
          if (per.le.acc) goto 20
!*
10     continue
!*
!** If fail to get required tolerence print warning.
!*
       write (nout,*) 'Ran out of iterations on bottom grid'
       write (nout,*) 'Increase itmax. itmax currently ',itmax
20     continue
!*
       dt = dthold
!*
!*
       return 
       end! Sroutine. 
!*
!**************** This is 72 columns wide *******************************
!*
       subroutine key(nlev,nmax,nl,ncx,narr)
!*
!* narr - length of grid array within xb
!* nlev - current level (labled from coursest(1) to finest grid)
!* nmax - maximum level used
!* nl   - starting point in arry
!* ncx  - current nx at this level
!*
       nl=1
!*
       do 10 l=1,nmax-nlev-1
          nx = 2**(nmax-l)
          nl = nl+(2*nx+1)**2
10     continue
!*
       ncx  = 2**nlev
       narr = (2*ncx+1)**2
!*
       return
       end ! key.
!*
!**************** This is 72 columns wide *******************************
!*










! *
       subroutine smoother(ffc,xb,xn,res,lxb,nx,dx)
       use sysdat_mod
       use times_mod
       use output_mod
       use ina_mod
! ***************************************************************
! *                                                             *
! * Smoother on all grids - Gauss-Seidel & SOR                  *
! *                                                             *
! ***************************************************************
! *
! *******This line is 72 columns wide*************************************
! *
	implicit double precision(a-h,o-z)
! *
        dimension xb(lxb),xn(lxb),ffc(lxb),res(lxb)
! *

       dt=0.25d0*dx*dx

! *
       call replace(xn,xb,lxb)
! *
       do 30 i = -nx,nx
       do 30 j = -nx,nx
! *
          call mfinder(nx,i,j,m)
          mip=m+1
          mim=m-1 
          mjp=m+2*nx+1
          mjm=m-2*nx-1 
! *
          if (i.eq.-nx.and.j.eq.-nx) then
             xn(m) =  0.5d0*(xn(mip)+xn(mjp)) - ffc(m)*dt
          elseif (i.eq.-nx.and.j.eq.nx) then
             xn(m) =  0.5d0*(xn(mip)+xn(mjm)) - ffc(m)*dt
          elseif (i.eq.-nx) then
             xn(m) = 0.5d0*xn(mip) + 0.25d0*(xn(mjp)+xn(mjm))&
                    -ffc(m)*dt
          elseif (i.eq.nx.and.j.eq.-nx) then
             xn(m) =  0.5d0*(xn(mim)+xn(mjp)) - ffc(m)*dt
          elseif (i.eq.nx.and.j.eq.nx) then
             xn(m) =  0.5d0*(xn(mim)+xn(mjm)) - ffc(m)*dt
          elseif (i.eq.nx) then
             xn(m) = 0.5d0*xn(mim) + 0.25d0*(xn(mjp)+xn(mjm))&
                    -ffc(m)*dt
          elseif (j.eq.-nx) then
             xn(m) = 0.5d0*xn(mjp) + 0.25d0*(xn(mip)+xn(mim))&
                    -ffc(m)*dt
          elseif (j.eq.nx) then
             xn(m) = 0.5d0*xn(mjm) + 0.25d0*(xn(mip)+xn(mim))&
                    -ffc(m)*dt
          else
             xn(m) = 0.25d0*( xn(mip)+xn(mim)+xn(mjp)+xn(mjm) )&
                    -ffc(m)*dt
          endif
! *
! *
30     continue
! *
! ** Calculate residuals
! *
! c      if (ina.eq.1.or.ina.eq.3) then
! *
! ** In interior
! *
          do 40 i=1-nx,nx-1
          do 40 j=1-nx,nx-1 
! *
             call mfinder(nx,i,j,m) 
             mip=m+1 
             mim=m-1  
             mjp=m+2*nx+1 
             mjm=m-2*nx-1  
! *
             xlhs =( xn(mip)-2.d0*xn(m)+xn(mim)&
                   +xn(mjp)-2.d0*xn(m)+xn(mjm) )/dx/dx
! *
             res(m) = ffc(m)-xlhs
40       continue
! *
! ** Boundary conditions
! *
       do 50 i= 1-nx,nx-1
          call mfinder(nx,i,-nx,m)
          mip=m+1
          mim=m-1 
          mjp=m+2*nx+1
! *
          xlhs =( xn(mip)-2.d0*xn(m)+xn(mim)&
                +2.d0*xn(mjp)-2.d0*xn(m) )/dx/dx
! *
          res(m) = ffc(m)-xlhs
! *
          call mfinder(nx,i, nx,m)
          mip=m+1
          mim=m-1
          mjm=m-2*nx-1
! *
          xlhs =( xn(mip)-2.d0*xn(m)+xn(mim)&
                +2.d0*xn(mjm)-2.d0*xn(m) )/dx/dx
! *
          res(m) = ffc(m)-xlhs
! *
          call mfinder(nx,-nx,i,m)
          mip=m+1
          mjp=m+2*nx+1
          mjm=m-2*nx-1
! *
          xlhs =( 2.d0*xn(mip)-2.d0*xn(m)&
                +xn(mjp)-2.d0*xn(m)+xn(mjm) )/dx/dx
! *
          res(m) = ffc(m)-xlhs
! *
          call mfinder(nx, nx,i,m)
          mim=m-1
          mjp=m+2*nx+1
          mjm=m-2*nx-1
! *
          xlhs =( 2.d0*xn(mim)-2.d0*xn(m) &
          +xn(mjp)-2.d0*xn(m)+xn(mjm) )/dx/dx
! *
          res(m) = ffc(m)-xlhs
50    continue
! *
! **  Corners
! *
      call mfinder(nx,-nx,-nx,m)
      mip=m+1
      mjp=m+2*nx+1
      xlhs = ( 2.d0*(xn(mip)+xn(mjp)) - 4.d0*xn(m) )/dx/dx
      res(m) = ffc(m)-xlhs
! *
      call mfinder(nx,-nx, nx,m)
      mip=m+1
      mjm=m-2*nx-1
      xlhs = ( 2.d0*(xn(mip)+xn(mjm)) - 4.d0*xn(m) )/dx/dx
      res(m) = ffc(m)-xlhs
! *
      call mfinder(nx, nx,-nx,m)
      mim=m-1
      mjp=m+2*nx+1
      xlhs = ( 2.d0*(xn(mim)+xn(mjp)) - 4.d0*xn(m) )/dx/dx
      res(m) = ffc(m)-xlhs
! *
      call mfinder(nx, nx, nx,m)
      mim=m-1
      mjm=m-2*nx-1
      xlhs = ( 2.d0*(xn(mim)+xn(mjm)) - 4.d0*xn(m) )/dx/dx
      res(m) = ffc(m)-xlhs
! *
! c      endif
! *
       return
! *
       end ! smoother






! *
       subroutine fcal(s,ffc,res,xn,ln,nx,dx)
       use sysdat_mod
! *
       implicit double precision(a-h,o-z)

! *
       double precision res(ln),xn(ln),ffc(ln)
! *

! *
       do 10 i=-nx,nx
       do 10 j=-nx,nx
          call mfinder(nx,i,j,m)
          ffc(m)=s*res(m)        ! For linear problems
10     continue
! *
       return
! *
       end ! fcal
! *


! *
       subroutine restrict(xn,ln,xb,lb,nx)
! *
       implicit double precision(a-h,o-z)
! *
! *
       double precision xn(ln),xb(lb)
! *
! * xn - input array to be restricted (fine grid).
! * xb - Output array (course grid).
! * nx - restricted grid size (course grid).
! * ln - length of xn
! * lb - length of xb
! *
! c      print*,'Restriction'
! *
! ** Loop to restric interior grid.
! *
! ** Full Weighting
! *
       do 10 i = 1-nx,nx-1
       do 10 j = 1-nx,nx-1 
          call mfinder(  nx,i,j,m1)        ! course grid
          call mfinder(2*nx,2*i,2*j,m2)    ! fine grid
! *                                          ! i,j
          call weight(w,xn,ln,m2,2*nx)
          xb(m1) = w/16.d0
 10    continue
! *
! ** Boundaries (Weighted)
! *
       do 30 i = 1-nx,nx-1
          call mfinder(nx,i,-nx,m1)        ! course grid
          call mfinder(2*nx,2*i,-2*nx,m2)  ! fine grid
! *
          m2a=m2+(4*nx+1)
          w = 2.d0*( xn(m2)+xn(m2a) )&
           + xn(m2-1)+xn(m2+1)+xn(m2a-1)+xn(m2a+1)
! *
          xb(m1) = w/8.d0
! *
          call mfinder(nx,i, nx,m1)        ! course grid
          call mfinder(2*nx,2*i, 2*nx,m2)  ! fine grid 
! *                                          ! i,j 
          m2a=m2-(4*nx+1) 
          w = 2.d0*( xn(m2)+xn(m2a) )&
           + xn(m2-1)+xn(m2+1)+xn(m2a-1)+xn(m2a+1)
! *
          xb(m1) = w/8.d0
 30    continue
! *
       do 40 j = 1-nx,nx-1
          call mfinder(nx,-nx,j,m1)        ! course grid
          call mfinder(2*nx,-2*nx,2*j,m2)  ! fine grid 
! *                                          ! i,j 
          m2p=m2+(4*nx+1)  
          m2m=m2-(4*nx+1)   
          w = 2.d0*( xn(m2)+xn(m2+1) ) &
           + xn(m2p)+xn(m2p+1)+xn(m2m)+xn(m2m+1) 
! *
          xb(m1) = w/8.d0
! * 
          call mfinder(nx, nx,j,m1)        ! course grid
          call mfinder(2*nx, 2*nx,2*j,m2)  ! fine grid  
! *                                          ! i,j  
          m2p=m2+(4*nx+1)   
          m2m=m2-(4*nx+1)    
          w = 2.d0*( xn(m2)+xn(m2-1) )  &
           + xn(m2p)+xn(m2p-1)+xn(m2m)+xn(m2m-1) 
! *
          xb(m1) = w/8.d0 
! *
 40    continue
! *
! ** Corners
! *
       call mfinder(nx,-nx,-nx,m1)       ! course grid 
       call mfinder(2*nx,-2*nx,-2*nx,m2) ! fine grid  
! *                                        ! i,j  
       m2p=m2+(4*nx+1)
       w=xn(m2)+xn(m2p)+xn(m2+1)+xn(m2p+1)
! *
       xb(m1) = w/4.d0 
! *
       call mfinder(nx,-nx,nx,m1)        ! course grid  
       call mfinder(2*nx,-2*nx,2*nx,m2)  ! fine grid  
! *                                        ! i,j   
       m2p=m2-(4*nx+1) 
       w=xn(m2)+xn(m2p)+xn(m2+1)+xn(m2p+1) 
! *
       xb(m1) = w/4.d0  
! *
       call mfinder(nx,nx,-nx,m1)        ! course grid  
       call mfinder(2*nx,2*nx,-2*nx,m2)  ! fine grid  
! *                                        ! i,j   
       m2p=m2+(4*nx+1) 
       w=xn(m2)+xn(m2p)+xn(m2-1)+xn(m2p-1) 
! *
       xb(m1) = w/4.d0  
! *
       call mfinder(nx,nx,nx,m1)         ! course grid  
       call mfinder(2*nx,2*nx,2*nx,m2)   ! fine grid  
! *                                        ! i,j   
       m2p=m2-(4*nx+1) 
       w=xn(m2)+xn(m2p)+xn(m2-1)+xn(m2p-1) 
! *
       xb(m1) = w/4.d0  
! *
       return
! *
       end ! restrictor.
! *
! **************************************************************************
! *
       subroutine weight(w,xn,ln,m2,nx)
! *
       implicit double precision(a-h,o-z)
! *
       dimension xn(ln)
! *
       mhold = m2         ! i,j
! *
       w  = 4.d0*xn(m2)
! *
       m2 = m2-1          ! i-1,j
       w  = w+2.d0*xn(m2)
! *
       m2 = m2+2          ! i+1,j
       w  = w+2.d0*xn(m2) 
! *
       m2 = m2-(2*nx+1)   ! i+1,j-1
       w  = w+xn(m2)
! *
       m2 = m2-1          ! i,j-1
       w  = w+2.d0*xn(m2)
! *
       m2 = m2-1          ! i-1,j-1
       w  = w+xn(m2)
! *
       m2 = m2+2*(2*nx+1) ! i-1,j+1
       w  = w+xn(m2)
! *
       m2 = m2+1          ! i,j+1
       w  = w+2.d0*xn(m2)
! *
       m2 = m2+1          ! i+1,j+1
       w  = w+xn(m2)
! *
       m2 = mhold       
! *
       return
       end ! weight.
! *
! **************************************************************************
! *
       subroutine mfinder(nx,i,j,m)
 
       implicit double precision (a-h,o-z) 
! *
       m=(j+nx)*(2*nx+1)+i+nx+1
! *
       return
       end ! mfinder.
! *





! *
       subroutine prolongate (xn,ln,xb,lb,nx)
       use sysdat_mod
! *
       implicit double precision(a-h,o-z)
! *
! *
       double precision xn(ln), xb(lb)
! *
! * nx - size of course grid
! * xb - array on course grid
! * xn - array on fine grid 
! *
! c      print*,'Prolongation'
! *
! ** Fill in boundaries of grid and points common to both course & fine grid.
! *
       do 20 j=-nx,nx
       do 20 i=-nx,nx
          call mfinder(2*nx,2*i,2*j,m1)  ! fine grid 
          call mfinder(  nx,i,j,m2)      ! course grid
          xn(m1) = xb(m2)
20     continue
! *
! ** Fill in every other point in x direction.
! *
       do 40 j=-nx,nx
       do 40 i=-nx,nx-1
          call mfinder(2*nx,2*i+1,2*j,m)  ! fine grid   
          xn(m) = 0.5d0*(xn(m-1)+xn(m+1))
40     continue
! *
! ** Fill in every other point in y direction. 
! * 
       nx2 = 4*nx+1
       do 60 i=-2*nx,2*nx
       do 60 j=-nx,nx-1
          call mfinder(2*nx,i,2*j+1,m)  ! fine grid    
          xn(m) = 0.5d0*(xn(m-nx2)+xn(m+nx2))
60     continue  
! *
       return
! *
       end! prolongate.
! *

      


