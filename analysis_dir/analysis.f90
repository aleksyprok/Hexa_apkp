program analysis
!Takes the output from hexa and calculates the free magnetic energy, helicty and current density and writes it out to file.
!Then run analysis.pro to view the output
!or run quick_analysis.pro to quickly view the evolution of the helicity and free magnetic energy
use variables
implicit none


dir = '../hexaf90/run1' !directory where hexa _#####p files are stored
dirp= '../potenMPI/run1' !directory where potential fields are stored
root='run1' !root of filename of hexa files
rootp='poten' !root of potential files

open (unit=10,file=trim(dir)//'/params.dat')
read (10,*) nmajor, nmajor
 !read in number of frames
close(10)


! Read in the param1 file
  fname=dir(1:length(dir,20))//'/param1'
  print *,'Read grid parameters: ',fname
  open(unit=3,file=fname,form='formatted',status='old')
  read (3,*) num_hex_cells
  read (3,*) nx,ny,nz
  read(3,*) l0
  close(3)

 print *, 'Grid sizes:', nx, ny, nz

 delx = 6./real(nx)
 dely = 6./real(nx)
 delz = 6./real(nx)



allocate(ax(nx+1,ny+1,nz+1),ay(nx+1,ny+1,nz+1),az(nx+1,ny+1,nz+1))
allocate(bx(nx+1,ny+1,nz+1),by(nx+1,ny+1,nz+1),bz(nx+1,ny+1,nz+1))
allocate(axp(nx+1,ny+1,nz+1),ayp(nx+1,ny+1,nz+1),azp(nx+1,ny+1,nz+1))
allocate(bxp(nx+1,ny+1,nz+1),byp(nx+1,ny+1,nz+1),bzp(nx+1,ny+1,nz+1))
allocate(cx(nx+1,ny+1,nz+1),cy(nx+1,ny+1,nz+1),cz(nx+1,ny+1,nz+1))
allocate(base(nx+1,ny+1,nmajor))
! apkp - start
allocate(bbx0(    nx + 1, 0 : ny + 1, nmajor))
allocate(bby0(0 : nx + 1,     ny + 1, nmajor))
! apkp - end

allocate(helicity(nmajor))


 open (unit=85,file=dir(1:length(dir,20))//'/analysis.out') !file containing basic properties of simulation (can be read in and plotted by quick_analysis.pro



 do nt=0,nmajor-1

   ftype = 'p'
   call readin(dirp,rootp) !readin potential fields


   ftype = 'h'
   call readin(dir,root) !readin hexa fields

   call calculate !calculate properties and write to file

enddo

 !write out helicity to file
 write (fname,fmt='(a,"/",a)')  dir(1:length(dir,20)),'helicity'
 print *,'Writing helicity to '//fname
     open(unit=56,file=fname,form='unformatted')
     write (56) nmajor
     write (56) (helicity)
     close(56)



 !write out base (bz on lower boundry) to file
 write (fname,fmt='(a,"/",a)')  dir(1:length(dir,20)),'base'
 open(unit=56,file=fname,form='unformatted')
     write (56) nmajor
     write (56) base
     close(56)

  ! apkp - start
  ! write out bbx0, bby0 on lower boundary to file
  write (fname,fmt='(a,"/",a)')  dir(1:length(dir,20)), 'bbx0_bby0'
  open(unit=56,file=fname,form='unformatted')
  write (56) bbx0
  write (56) bby0
  close(56)

  deallocate(bbx0, bby0)
  ! apkp - end

 end program
