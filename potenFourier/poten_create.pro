pro poten_create

openr,10,'params.dat'

siz=0
nt=0
alpha=0
dir='a'
periodic=0
open=0

readf,10,siz,nt,alpha,dir
readf,10,periodic,open
close,10

if (periodic eq 0) and (open eq 0) then begin
  print, 'Computational box has closed top and side boundaries.'
  print, 'Please use the codes in the potenMPI/ directory.'
  stop
endif

dir=strtrim(dir,2)

print,nt
print,dir
print,periodic
print,open


nx=2*2^siz
ny=nx
nz=nx

print,'nx=',nx
print,'ny=',ny
print,'nt=',nt

bz=dblarr(nx,ny,nt)

openu,10,dir+'mag_data.dat',/f77_unformatted
  readu,10,bz
close,10


imb=fltarr(nt)
for t=0,nt-1 do begin
  imb(t)=total(bz(*,*,t))/nx/ny
  bz(*,*,t)=bz(*,*,t)-imb(t)
endfor


window,0,xs=nx,ys=ny

delx = 6./nx
dely = delx
delz = delx

filename='poten'

if periodic eq 1 then begin
  xbc = 1       ; 0 = closed, 1 = periodic
  ybc = 1       ; 0 = closed, 1 = periodic
endif else begin
  xbc = 0       ; 0 = closed, 1 = periodic
  ybc = 0       ; 0 = closed, 1 = periodic
endelse

if open eq 1 then begin
  topbc = 1     ; 0 = closed, 1 = open
endif else begin
  topbc = 0
endelse


print,'x-periodic y-periodic open'
print,xbc,ybc,topbc

;***************************************



for t=0,nt-1 do begin ;loop over every magnetogram

    bz0 = bz(*,*,t)     ; select frame

    ;-----------------------------------------------------------
    ; Compute potential field:
    ;-----------------------------------------------------------

    if xbc eq 0 then nx2=2*nx else nx2=nx         ; 0 = closed (reflective), 1 = periodic
    if ybc eq 0 then ny2=2*ny else ny2=ny

    bzd=fltarr(nx2,ny2)

    kx=findgen(nx2)/nx2
    kx(nx2/2+1:nx2-1)=kx(nx2/2+1:nx2-1)-1.0
    ky=findgen(ny2)/ny2
    ky(ny2/2+1:ny2-1)=ky(ny2/2+1:ny2-1)-1.0

    qq=4./delx^2*((sin(!pi*kx)^2)#(fltarr(ny2)+1.0))   $
      +4./dely^2*((fltarr(nx2)+1.0)#(sin(!pi*ky)^2))
    qq(0,0)=1.0

    aax=fltarr(nx  ,ny+1,nz+1)
    aay=fltarr(nx+1,ny  ,nz+1)
    aaz=fltarr(nx+1,ny+1,nz  )

    bzd(0:nx-1, 0:ny -1)=bz0
    if xbc eq 0 then bzd(nx:nx2-1,0:ny-1)=reverse(bz0,1)    ; Closed in x (mirror in x)
    if ybc eq 0 then bzd(0:nx-1,ny:ny2-1)=reverse(bz0,2)    ; Closed in y (mirror in y)
    if ybc eq 0 and xbc eq 0 then bzd(nx:nx2-1,*)=reverse(bzd(0:nx-1,*),1) ; make sure mirrored in all 4 quadrants

    ht=fft(bzd,-1)/qq
    q=sqrt(double(qq))
    q(0,0)=0.d0

    ztop=double(nz*delz)
    denom=1.d0-exp(-((2.d0*ztop*q)<20.))
    denom(0,0)=1.d0

     
      for k=0,nz do begin
      
	z=double(k*delz)

    ;  If top boundary closed, Bz must decrease to zero at top:

	if (topbc eq 0) then $ 
	    fact=float((exp(-((q*z)<20.))-exp(-((q*(2.d0*ztop-z))<20.)))/denom) $     ; Closed top
	else fact=float(exp(-((q*z)<20.)))                                               ; Open top
	fact(0,0)=0.

	h=float(fft(fact*ht,1))
	aax(0:nx-1,1:ny-1,k)= (h(0:nx-1,1:ny-1)-h(0:nx-1,0:ny-2))/dely
	if (ybc eq 1) then begin
	  aax(0:nx-1,     0,k)= (h(0:nx-1,     0)-h(0:nx-1,  ny-1))/dely
	  aax(0:nx-1,    ny,k)=aax(0:nx-1,0,k)
	endif

	aay(1:nx-1,0:ny-1,k)=-(h(1:nx-1,0:ny-1)-h(0:nx-2,0:ny-1))/delx
	if (xbc eq 1) then begin
	  aay(     0,0:ny-1,k)=-(h(     0,0:ny-1)-h(  nx-1,0:ny-1))/delx
	  aay(    nx,0:ny-1,k)=aay(0,0:ny-1,k)
	endif
	
	for j=0,ny-1 do begin
	  aay(*,j,k)=aay(*,j,k)+imb(0)*findgen(nx+1)*delx
	endfor
      endfor


    ;-----------------------------------------------------------
    ; Save to dir/filename_#####p
    ;-----------------------------------------------------------

    fullname=string(filename,t,format='(a,"_",i5.5,"p")')
    file=dir+fullname
    print,'Save '+file
    openw,1,file,/f77_unform
    opt=1L
    ;print,string(opt,format='("File type: opt=",i2)')
    writeu,1,opt
    writeu,1,aax
    writeu,1,aay
    writeu,1,aaz
    close,1
    
    bz1=(aay(1:nx,*)-aay(0:nx-1,*))/delx - (aax(*,1:ny)-aax(*,0:ny-1))/dely

    tv,bytscl(bz1,-100,100)
    

endfor


print,'Done!'


end