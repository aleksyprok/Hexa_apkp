module var_global

  implicit none
  include "mpif.h"
  !use MPI
  real, parameter :: pi=3.141592654
  real, parameter :: frc_coef=3000. !Frictional coefficient (km^2/s)

!  Grid parameters:

  real    :: delx,dely,delz

!  Number of cells (not corners):

  integer :: nxglobal,nyglobal,nzglobal
  integer :: nx,ny,nz

!  Vector potentials:

  !integer, parameter :: ntmax=401
  integer :: ntmax
  real, dimension(:,:,:), allocatable :: aax
  real, dimension(:,:,:), allocatable :: aay
  real, dimension(:,:,:), allocatable :: aaz
  real, dimension(:,:,:), allocatable :: aax0
  real, dimension(:,:,:), allocatable :: aay0
  real, dimension(:,:), allocatable :: daax
  real, dimension(:,:), allocatable :: daay

  real, dimension(:,:,:), allocatable :: div,bx,by,bz,bb,bbm
  real, dimension(:,:,:), allocatable :: cx,cy,cz,ch,vx,vy,vz
  real, dimension(:,:,:), allocatable :: bbx
  real, dimension(:,:,:), allocatable :: bby
  real, dimension(:,:,:), allocatable :: bbz

  real, dimension(:,:,:), allocatable :: ccx
  real, dimension(:,:,:), allocatable :: ccy
  real, dimension(:,:,:), allocatable :: ccz

!  Model parameters:

  integer :: vsetup,nrelax,nmajor,nminor,nrep
  real :: bzcrit,dtminor,dtime,etad,etai,eta4

! Subdirectory name.


    character*30 :: dir
    character*20 :: root

    character*50 :: filename

! Input variable for run


  integer :: nstrt,nend,mode
  real :: etaia,eta4a
  integer :: nrun,opt
  integer :: periodic, open

! variable associated with MPI

integer :: mpisize, ierr, rank,comm,left,right,up,down,nextrank,rankstart,rankend
integer, DIMENSION(MPI_STATUS_SIZE) :: stat
integer, parameter :: mpidir=2
integer, dimension(mpidir) :: dims,nproc,coords
logical, dimension(mpidir) :: periods
integer, parameter :: tag=100,tag1=101,tag2=102

! variables associated with length and time

real :: length_cm !length of one Hexa length unit in cm
real :: time_s !time of one hexa time unit in seconds
real :: dt,basedt !timestep and default timestep (hexa units)
real :: timestep_s !default timestep in seconds
real :: time !time in hexa units


! variables associated with flux imbalance corrections

real, dimension(:), allocatable :: corr,fimb
real :: imb,totimb
real, dimension(:,:,:), allocatable :: ay_corr

  contains

!-----------------------------------------------------------------
!  Function "length" determines non-blank length of a string:
!-----------------------------------------------------------------
    function length(char,max)
      integer :: length,max
      character*(*) char

      length=max+1
 8    length=length-1
      if(length.eq.0) return
      if(char(length:length).eq.' ') goto 8

    end function length

!-----------------------------------------------------------------
!  Get parameter value:
!-----------------------------------------------------------------

   function get_value(unit,char)

    integer :: unit,i,n1,n2
    real :: get_value
    character*(*) char
    character*80 text
    character*20 name

    read (unit,*) text
!   print *,text

    n1=1
    do i=1,80
       if(text(i:i).eq.' ') n1=i+1
       if(text(i:i).eq.'=') then
          n2=i-1
          if(n2.lt.n1) then
             print *,'Error: No name before equal sign'
             goto 9
          endif
          goto 8
       endif
    enddo
    print *,'Error: Missing equal sign'
    goto 9
8   continue

!   print *,'n1=',n1,' n2=',n2

    name=text(n1:n2)
    if(name(1:length(name,20)).ne.char) then
       print *,'Error: Incorrect parameter name'
       print *,'Expected ',char,', found ',name
       goto 9
    endif
    i=n2+2
    if(length(text(i:80),80-i+1).eq.0) then
       print *,'Error: Missing parameter value'
       goto 9
    endif
    read (text(i:80),*) get_value
    return

9   print *,'*** Error reading parameter file:'
    print *,text
    stop 'Error'

  end function get_value

end module var_global
