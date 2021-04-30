module sysdat_mod
!Replaces the 'sysdat.dat' file
integer :: nmax, na, nb, nm
end



module avalues_mod
double precision, allocatable :: ax(:,:), ay(:,:)
end



module dunk_mod
integer :: dunk_count,nt
end



module xn_mod
double precision, allocatable :: xn(:)
end



module xb_mod
double precision, allocatable :: xb(:)
end



module times_mod
double precision :: dt, time
end



module output_mod
integer :: nout
end



module rhs1_mod
double precision, allocatable :: ff(:)
end



module ina_mod
integer :: ina
end



module wu_mod
double precision :: wu
end



module res_mod
double precision, allocatable :: res(:), resc(:)
end



module acc_mod
double precision :: acc
end



module const_mod
double precision :: dx, pi
integer :: nx
end



module xcb_mod
double precision, allocatable :: xcb(:), xcn(:), xct(:)
end



module rhs_mod
double precision, allocatable :: ffc(:)
end



module s_mod
double precision, allocatable :: s(:)
end



module cor_mag_mod
double precision, allocatable :: cor_mag(:,:,:)
end



module cdata
!The replacement of CDATA.IN in the old file
integer :: nx0
double precision ::alpha,H0
double precision, parameter :: OM=1.5
double precision, parameter :: Delta=1.d-8
integer, parameter :: nlev=7
double precision, parameter ::wmax=302.0
integer, parameter :: maxits=30
integer, parameter :: istart=0
end


module fnames
character (len=20) :: dir
end








subroutine array_alloc
!allocates arrays for the evolve routines
use sysdat_mod
use rhs1_mod
use s_mod
use rhs_mod
use xcb_mod
use xb_mod
use xn_mod
use res_mod
use dunk_mod
use cor_mag_mod

allocate(ff(nm))
allocate(s(nmax))
allocate(ffc(nm/2))
allocate(xcb(nm/2),xcn(nm/2),xct(nm/2))
allocate(xb(nm),xn(nm))
allocate(res(nm),resc(nm/2))

allocate(cor_mag(2*na,2*na,nt))
end



subroutine array_dealloc
!deallocates arrays allocated in array_alloc
use sysdat_mod
use rhs1_mod
use s_mod
use rhs_mod
use xcb_mod
use xb_mod
use xn_mod
use res_mod

deallocate(ff)
deallocate(s)
deallocate(ffc)
deallocate(xcb,xcn,xct)
deallocate(xb,xn)
deallocate(res,resc)

end





subroutine sysdat
!Reads in gridsize (2*2**nmax) and the number of frames (nt)
!replaces sysdat.dat in the old code
use sysdat_mod
use dunk_mod
use cdata
use fnames

character (len=60) :: fmt
integer :: periodic,open

open (unit=10,file='params.dat')
read(10,*) nmax, nt, alpha,dir
read(10,*) periodic,open

if (periodic .eq. 1 .or. open .eq. 1) then
  print*, 'Incompatible boundary conditions (open or periodic selected)'
  print*, 'Please use the codes in the potenFourier/ directory'
  stop
endif

na=2**nmax
nm=(2*na+1)**2 
nb=na

nx0=2*na+1
h0=2./dble(nx0-1)



print*,'Working directory: ',trim(dir)
fmt="(' Gridsize= ',i3,'^3')"
write(*,fmt) nx0-1
fmt="(' Number of frames=',i4.1)"
write(*,fmt) nt
fmt="(' Alpha=',F5.3)"
!print*,'ALPHA=',alpha
write(*,fmt) alpha



close(10)

end




subroutine readin
!reads in cor_mag from 'mag_data'
use cor_mag_mod
use sysdat_mod
use dunk_mod
use fnames



open (unit=10,file=trim(dir)//'/mag_data.dat',form='unformatted')
read(10),cor_mag
close(10)

end






subroutine writeout_evolve
!writes Ax and Ay out to a hexa '_evolve' file
use sysdat_mod
use avalues_mod

real :: aax(2*na,2*na+1),aay(2*na+1,2*na)

!aax(:,:) =0.5*real(ax(-na:na-1,:)+ax(1-na:na,:))*3.
!aay(:,:) = 0.5*real(ay(:,-na:na-1) + ay(:,1-na:na))*3.

aax(:,:) = 1.5*real(ax(1:2*na,:)+ax(2:2*na+1,:))
aay(:,:) = 1.5*real(ay(:,1:2*na)+ay(:,2:2*na+1))
write(11) aax
write(11) aay

end


subroutine writeout_whole2(ax,ay,az)
use cdata
use fnames
use dunk_mod
double precision :: ax(nx0,nx0,nx0),ay(nx0,nx0,nx0),az(nx0,nx0,nx0)
integer :: opt
real :: aax(nx0-1,nx0,nx0), aay(nx0,nx0-1,nx0), aaz(nx0,nx0,nx0-1)
real :: dx
character (len=50) :: filename

!  OPEN(UNIT=10,FILE='AXOUT',FORM='UNFORMATTED')
!       WRITE(10) Ax
!       CLOSE(10)
!       OPEN(UNIT=10,FILE='AYOUT',FORM='UNFORMATTED')
!       WRITE(10) Ay
!       CLOSE(10)
!       OPEN(UNIT=10,FILE='AZOUT',FORM='UNFORMATTED')
!       WRITE(10) Az
!       CLOSE(10)
!       

!average As to cell corners. Multiply by 3 to get in Hexa units
aax(:,:,:) =3.*(ax(1:nx0-1,:,:)+ax(2:nx0+1,:,:))/2.
aay(:,:,:) =3.*(ay(:,1:nx0-1,:)+ay(:,2:nx0+1,:))/2.
aaz(:,:,:) =3.*(az(:,:,1:nx0-1)+az(:,:,2:nx0+1))/2.

opt=1

write (filename,fmt='("poten_",i5.5,"p")') dunk_count-1
      
      open(unit=10,file=trim(dir)//'/'//filename,form='unformatted')

!open (unit=10,file=trim(dir)//'/hexa_files/'//trim(dir)//'_00000p',form='unformatted')
write(10) opt
write(10) aax
write(10) aay
write(10) aaz
close(10)

print*,'File ',trim(trim(dir)//'/'//filename),' written to file.'

end



     SUBROUTINE INITIALAS(NX,NY,NZ,H,AX2,AY2)
     !takes Ax and Ay from poislin and feeds them into the 3D Ax and Ay arrays for f3cray
      
      use avalues_mod
      implicit none
      INTEGER NX,NY,NZ
      DOUBLE PRECISION H,AX2(NX,NY,NZ),AY2(NX,NY,NZ)
      INTEGER I,J

      ax2(:,:,1)=ax
      ay2(:,:,1)=ay
       
     ! deallocate(ax,ay)
      
      END





module duncan
double precision, allocatable :: ax(:), ay(:), az(:)
end



subroutine transmit
use mpi
use sysdat_mod
use dunk_mod
use cdata
use fnames

call MPI_Bcast(nmax,1,MPI_INT,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(nt,1,MPI_INT,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
call MPI_Bcast(dir,20,MPI_CHAR,0,MPI_COMM_WORLD,ierror)

na=2**nmax
nm=(2*na+1)**2 
nb=na

nx0=2*na+1
h0=2./dble(nx0-1)

end