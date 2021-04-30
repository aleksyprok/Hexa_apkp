       subroutine bzfinder()
!Takes magnetogram data (in cor_mag) and feeds it into the 1D array ff for
!use in the calculations
       use dunk_mod
       use sysdat_mod
       use const_mod
       use rhs1_mod
       use cor_mag_mod
       use output_mod

       implicit double precision(a-h,o-z)

!feed cor_mag into ff
       do 11 j = 1-nx,nx-1
       do 11 i = 1-nx,nx-1
          test = 0
          call mfinder(nx,i,j,m)
             ff(m) = -cor_mag(i+nx,j+nx,dunk_count)
11     continue
      
! The below steps ensure flux balancing


!Determine flux through base
       flux = 0.d0
       do 150 j=1-nx,nx-1
       do 150 i=1-nx,nx-1
          call mfinder(nx,i,j,m)
          flux = flux + ff(m)
 150   continue



!determine # of elements where |B|>1G in magnetogram
       icounter=0
       do 151 j=1-nx,nx-1
         do 151 i=1-nx,nx-1
           call mfinder(nx,i,j,m)
           if (dabs(ff(m)) .gt. 1.d0) icounter=icounter+1
 151   continue
       flux = flux/dble(icounter) !flux correction per pixel needed
       if (abs(flux) .gt. 1.) then
         print*,'Flux correction changes pixel signs!',dunk_count,flux
       endif
       
       
! c subtract flux correction from pixels     
       do 152 j=1-nx,nx-1
         do 152 i=1-nx,nx-1
           call mfinder(nx,i,j,m)
           if (dabs(ff(m)) .gt. 1.d0) ff(m) = ff(m) - flux
 152   continue
        

       return

       end 

