module variables
implicit none

integer :: nmajor


double precision, allocatable, dimension(:,:,:) :: ax, ay, az, bx, by, bz
double precision, allocatable, dimension(:,:,:) :: axp, ayp, azp, bxp, byp, bzp
double precision, allocatable :: helicity(:)
double precision :: L0
double precision, allocatable, dimension(:,:,:) :: cx, cy, cz
double precision, allocatable, dimension(:,:,:) :: base
! apkp - start
double precision, allocatable, dimension(:,:,:) :: bbx0, bbz0
! apkp - end

 character (len=40) :: dir, root,dirp,rootp
 character (len=40) :: fname
 character (len=1) :: ftype
 integer :: opt, i, j, k, nx, ny, nz, num_hex_cells, nt
 real :: delx, dely, delz

contains
      function length(char,max)
      integer :: length,max
      character*(*) char

      length=max+1
 8    length=length-1
      if(length.eq.0) return
      if(char(length:length).eq.' ') goto 8

    end function length
end module
