subroutine calculate
use variables
implicit none

double precision,dimension(nx+1,ny+1,nz+1) :: efree, bp2, b2
double precision :: hel(nx+1,ny+1,nz+1)

! Free Energy Calculation

b2 = (bx*bx + by*by + bz*bz)
bp2 = (bxp*bxp + byp*byp + bzp*bzp)

efree=b2-bp2
efree = efree/8.d0/3.1415926535d0*(L0*delx)**3
b2=b2/8.d0/3.1415926535d0*(L0*delx)**3
bp2=bp2/8.d0/3.1415926535d0*(L0*delx)**3

 cx = cx*cx+cy*cy+cz*cz

write (fname,fmt='(a,"/",a,"_",i5.5,"p")')  dir(1:length(dir,20)),'fenergy',nt
     print *,'Writing free energy to '//fname
     open(unit=56,file=fname,form='unformatted')
     write (56) nx,ny,nz
     write (56) (efree)
     write(56) b2
     write(56) cx
 
     close(56)

!A.B - Ap.Bp   (berger + field 1984)  
hel= ax*bx+ay*by+az*bz 
hel=hel-(axp*bxp+ayp*byp+azp*bzp)

!Ap.B - A.Bp (Finn and Antonsen 1985)
hel=hel+axp*bx+ayp*by+azp*bz
hel=hel-(ax*bxp+ay*byp+az*bzp)

helicity(nt+1) = sum(hel)




helicity(nt+1) = helicity(nt+1)*(delx**3)*L0**4 !NB the extra L0 comes from the fact that dimensionally A=B*L0
print*,'helicity', helicity(nt+1)

write(85,*) nt,helicity(nt+1),sum(b2),sum(bp2),sum(efree)
call flush(85)
end
