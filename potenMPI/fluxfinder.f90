       subroutine bzfinder()
       
       use dunk_mod
       use sysdat_mod
       use const_mod
       use rhs1_mod
       use cor_mag_mod
        use output_mod
!*
       implicit double precision(a-h,o-z)

       do 11 j = 1-nx,nx-1
       do 11 i = 1-nx,nx-1
          test = 0
          call mfinder(nx,i,j,m)
             ff(m) = -cor_mag(i+nx,j+nx,dunk_count)
11     continue
      


       flux = 0.d0
       do 150 j=1-nx,nx-1
       do 150 i=1-nx,nx-1
          call mfinder(nx,i,j,m)
          flux = flux + ff(m)
 150   continue
     !  print*,'Flux before correction',flux



       icounter=0
       do 151 j=1-nx,nx-1
         do 151 i=1-nx,nx-1
           call mfinder(nx,i,j,m)
           if (dabs(ff(m)) .gt. 1.d0) icounter=icounter+1
 151   continue
       flux = flux/dble(icounter)
       
       
! c subtrct flux (Gordon)       
       do 152 j=1-nx,nx-1
         do 152 i=1-nx,nx-1
           call mfinder(nx,i,j,m)
           if (dabs(ff(m)) .gt. 1.d0) ff(m) = ff(m) - flux
 152   continue
        

       flux = 0.d0
       do 155 j=1-nx,nx-1
       do 155 i=1-nx,nx-1
          call mfinder(nx,i,j,m)
          flux = flux + ff(m)
 155   continue
      ! print*,'Flux after correction',flux
       
! 
! **
! ** Sides
       do 160 i = 1-nx,nx-1
          call mfinder(nx,i,-nx,m)
          flux = flux + 0.5d0*ff(m)
          call mfinder(nx,i, nx,m)
          flux = flux + 0.5d0*ff(m)
          call mfinder(nx,-nx,i,m)
          flux = flux + 0.5d0*ff(m)
          call mfinder(nx, nx,i,m)
          flux = flux + 0.5d0*ff(m)
 160   continue
! *
! ** Corners
! *
       call mfinder(nx,-nx,-nx,m)
       flux = flux + 0.25d0*ff(m)
       call mfinder(nx,-nx, nx,m)
       flux = flux + 0.25d0*ff(m)
       call mfinder(nx, nx,-nx,m)
       flux = flux + 0.25d0*ff(m)
       call mfinder(nx, nx, nx,m)
       flux = flux + 0.25d0*ff(m)
! *
! *       write(nout,*) ' '
  !     write(nout,997) flux
! *

!997    format (2x,'Total flux through base: ',d12.6)
!998    format (8x,'x_0',11x,'y_0',11x,'re',11x,'B_0')
!999    format (4(2x,f12.6),i3)
 

     
! *
       return
! *
       end 

