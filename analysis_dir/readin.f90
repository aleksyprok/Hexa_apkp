subroutine readin(dir2,root2)
use variables
implicit none

!real, allocatable :: aax(:,:,:), aay(:,:,:), aaz(:,:,:)

real :: aax1(nx,ny+1,nz+1),aay1(nx+1, ny, nz+1),aaz1(nx+1,ny+1,nz)
double precision :: aax(nx,ny+1,nz+1),aay(nx+1, ny, nz+1),aaz(nx+1,ny+1,nz)
double precision :: bbx(nx+1,0:ny+1,0:nz+1),bby(0:nx+1,ny+1,0:nz+1),bbz(0:nx+1,0:ny+1,nz+1)
double precision :: ccx(nx+2,ny+1,nz+1),ccy(nx+1,ny+2,nz+1),ccz(nx+1,ny+1,nz+2)
character (len=40) :: dir2,root2

 !allocate(aax(nx,ny+1,nz+1),aay(nx+1, ny, nz+1),aaz(nx+1,ny+1,nz))

 if (ftype .ne. 'p' .and. ftype .ne. 'h') then
   print *, "Incorrect Read Parameter"
   stop
 endif


 write (fname,fmt='(a,"/",a,"_",i5.5,"p")')  &     ! *** Mac
          dir2(1:length(dir2,20)),root2(1:length(root2,20)),nt
          print*, fname
     !print *,'Reading 3D model from '//fname
     open(unit=56,file=fname,form='unformatted',status='old')
     read (56) opt
     !print *,'File type: opt=',opt
     if(.not.(opt.eq.1)) stop 'Invalid option'
     read (56) aax1
     read (56) aay1
     read(56) aaz1
     close(56)

 aax = dble(aax1)
aay = dble(aay1)
aaz = dble(aaz1)



!calculate ax,ay,az at cell corner:

if (ftype .eq. 'h') then
!print *, 'Calculating grid corner values - non-potential'


  ax(2:nx,:,:)=0.5*(aax(1:nx-1,:,:)+aax(2:nx,:,:))
  ay(:,2:ny,:)=0.5*(aay(:,1:ny-1,:)+aay(:,2:ny,:))
  az(:,:,2:nz)=0.5*(aaz(:,:,1:nz-1)+aaz(:,:,2:nz))



else if (ftype .eq. 'p') then
!print *, 'Calculating grid corner values - potential'
  axp(2:nx,:,:)=0.5*(aax(1:nx-1,:,:)+aax(2:nx,:,:))
  ayp(:,2:ny,:)=0.5*(aay(:,1:ny-1,:)+aay(:,2:ny,:))
  azp(:,:,2:nz)=0.5*(aaz(:,:,1:nz-1)+aaz(:,:,2:nz))


else
  print *, 'Incorrect type'
  stop
endif


!print *, 'calculating b on faces'

bbz(1:nx,1:ny,1:nz+1)=  &
       (aay(2:nx+1,1:ny,1:nz+1)-aay(1:nx,1:ny,1:nz+1))/delx  &
      -(aax(1:nx,2:ny+1,1:nz+1)-aax(1:nx,1:ny,1:nz+1))/dely


bbx(1:nx+1,1:ny,1:nz)=  &
       (aaz(1:nx+1,2:ny+1,1:nz)-aaz(1:nx+1,1:ny,1:nz))/dely  &
      -(aay(1:nx+1,1:ny,2:nz+1)-aay(1:nx+1,1:ny,1:nz))/delz


bby(1:nx,1:ny+1,1:nz)=  &
       (aax(1:nx,1:ny+1,2:nz+1)-aax(1:nx,1:ny+1,1:nz))/delz  &
      -(aaz(2:nx+1,1:ny+1,1:nz)-aaz(1:nx,1:ny+1,1:nz))/delx

!print *, 'max,bz', maxval(bbz(:,:,1))

!Boundary conditions:
!#################################
bbz(0,:,:) = bbz(1,:,:)
bbz(nx+1,:,:) = bbz(nx,:,:)

bbz(:,0,:) = bbz(:,1,:)
bbz(:,ny+1,:) = bbz(:,ny,:)
!##################################
bbx(:,0,:) = bbx(:,1,:)
bbx(:,ny+1,:) = bbx(:,ny,:)

bbx(:,:,0) = bbx(:,:,1)

 bbx(:,:,   0)=bbx(:,:,   1)  &
         -delz/delx*(bbz(1:nx+1,0:ny+1,1)-bbz(0:nx+0,0:ny+1,1))

bbx(:,:,nz+1) = bbx(:,:,nz)
!##################################
bby(0,:,:) = bby(1,:,:)

bby(:,:,   0)=bby(:,:,   1)  &
         -delz/dely*(bbz(0:nx+1,1:ny+1,1)-bbz(0:nx+1,0:ny+0,1))

bby(nx+1,:,:) = bby(nx,:,:)

bby(:,:,0) = bby(:,:,1)
bby(:,:,nz+1) = bby(:,:,nz)
!##################################


if (ftype .eq. 'h') then


bz=(bbz(0:nx,0:ny,:)+bbz(1:nx+1,0:ny,:)+bbz(0:nx,1:ny+1,:)+bbz(1:nx+1,1:ny+1,:))/4.

bx=(bbx(:,0:ny,0:nz)+bbx(:,1:ny+1,0:nz)+bbx(:,0:ny,1:nz+1)+bbx(:,1:ny+1,1:nz+1))/4.

by=(bby(0:nx,:,0:nz)+bby(1:nx+1,:,0:nz)+bby(0:nx,:,1:nz+1)+bby(1:nx+1,:,1:nz+1))/4.


!Calculate current on edges
    ccz(1:nx+1,1:ny+1,1:nz+2)=  &
         (bby(1:nx+1,1:ny+1,0:nz+1)-bby(0:nx,1:ny+1,0:nz+1))/delx  &
        -(bbx(1:nx+1,1:ny+1,0:nz+1)-bbx(1:nx+1,0:ny,0:nz+1))/dely

    ccx(1:nx+2,1:ny+1,1:nz+1)=  &
         (bbz(0:nx+1,1:ny+1,1:nz+1)-bbz(0:nx+1,0:ny,1:nz+1))/dely  &
        -(bby(0:nx+1,1:ny+1,1:nz+1)-bby(0:nx+1,1:ny+1,0:nz))/delz

    ccy(1:nx+1,1:ny+2,1:nz+1)=  &
         (bbx(1:nx+1,0:ny+1,1:nz+1)-bbx(1:nx+1,0:ny+1,0:nz))/delz  &
        -(bbz(1:nx+1,0:ny+1,1:nz+1)-bbz(0:nx,0:ny+1,1:nz+1))/delx



    cx=0.5*(ccx(1:nx+1,:,:)+ccx(2:nx+2,:,:))
    cy=0.5*(ccy(:,1:ny+1,:)+ccy(:,2:ny+2,:))
    cz=0.5*(ccz(:,:,1:nz+1)+ccz(:,:,2:nz+2))






base(:,:,nt+1) = bz(:,:,1)
! apkp - start
bbx0(:, :, nt + 1) = 0.5 * (bbx(:, :, 0) + bbx(:, :, 1))
bby0(:, :, nt + 1) = 0.5 * (bby(:, :, 0) + bby(:, :, 1))
! apkp - end

else if (ftype .eq. 'p') then

bzp=(bbz(0:nx,0:ny,:)+bbz(1:nx+1,0:ny,:)+bbz(0:nx,1:ny+1,:)+bbz(1:nx+1,1:ny+1,:))/4.

bxp=(bbx(:,0:ny,0:nz)+bbx(:,1:ny+1,0:nz)+bbx(:,0:ny,1:nz+1)+bbx(:,1:ny+1,1:nz+1))/4.

byp=(bby(0:nx,:,0:nz)+bby(1:nx+1,:,0:nz)+bby(0:nx,:,1:nz+1)+bby(1:nx+1,:,1:nz+1))/4.


else
  print *, 'ftype error'
  stop
endif







end subroutine
