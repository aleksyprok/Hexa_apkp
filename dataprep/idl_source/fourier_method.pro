pro fourier_method, bz_mag,dir,filename,periodic,open

;
;  Compute potential field and evolve file from a magnetogram series
;  called bz_mag, of size nx,ny,nt
;



siz = size(bz_mag)         ; determine magnetogram dimensions
nx = siz(1)
ny = siz(2)
nt = siz(3)


;bz_mag=congrid(bz_mag,256,256,nt)

; bz_mag(0,*,*)=0
; bz_mag(255,*,*)=0
; bz_mag(*,255,*)=0
; bz_mag(*,0,*)=0

;nx=256
;ny=256




print,nx,ny,nt


imb=fltarr(nt)

;calculate flux imbalance and correct
; this is needed to solve the equations
; (flux is 're-unbalanced' later)
for t=0,nt-1 do begin
  imb(t)=total(bz_mag(*,*,t))/nx/ny
  ;imb(t)=1.0
  bz_mag(*,*,t)=bz_mag(*,*,t)-imb(t)
endfor
;imb=imb/2



;setup gridsize
nz   = nx
delx = 6./nx
dely = delx
delz = delx
; dtime=1./(5000)
; eta4=0.001*(delx^4)/dtime


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

bz0 = bz_mag(*,*,0)     ; Initial frame



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


;shade_surf,bzd,charsize=2,xtitle='x',ytitle='y'
;stop

ht=fft(bzd,-1)/qq
q=sqrt(double(qq))
q(0,0)=0.d0

ztop=double(nz*delz)
denom=1.d0-exp(-((2.d0*ztop*q)<20.))
denom(0,0)=1.d0

  print,'Computing potential field:'
  for k=0,nz do begin
   ; if (k mod 5) eq 0 then print,string(k,format='("height =",i3)')
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
      aay(*,j,k)=aay(*,j,k)+imb(0)*(findgen(nx+1))*delx
    endfor
    
  endfor

  print,'done!'
;-----------------------------------------------------------
; Save to dir/filename_00000p
;-----------------------------------------------------------

fullname=string(filename,0,format='(a,"_",i5.5,"p")')
file=dir+'hexa_files/'+fullname
print,'Save '+file
openw,1,file,/f77_unform
opt=1L
print,string(opt,format='("File type: opt=",i2)')
writeu,1,opt
writeu,1,aax
writeu,1,aay
writeu,1,aaz
close,1





;stop

print,''
print,''

;-----------------------------------------------------------
; Compute Ax and Ay for photospheric boundary evolution:
;-----------------------------------------------------------

print,'Computing lower boundary:'

aax0=fltarr(nx  ,ny+1)
aay0=fltarr(nx+1,ny  )
file=dir+'hexa_files/'+filename+'_evolve'

openw,1,file,/f77_unform
opt=1L
writeu,1,opt
print,'Opening '+file

axg=fltarr(nx,ny+1,nt)
ayg=fltarr(nx+1,ny,nt)
window,0,xs=nx,ys=ny






for t=0,nt-1 do begin
aax0(*,*)=0
aay0(*,*)=0


  bzd(0:nx-1, 0:ny -1)=bz_mag(*,*,t)
  if xbc eq 0 then bzd(nx:nx2-1,0:ny-1)=reverse(bzd(0:nx-1,0:ny-1),1)    ; Closed in x
  if ybc eq 0 then bzd(0:nx-1,ny:ny2-1)=reverse(bzd(0:nx-1,0:ny-1),2)    ; Closed in y
  if (ybc eq 0) and (xbc eq 0) then bzd(nx:nx2-1,*)=reverse(bzd(0:nx-1,*),1) ; make sure mirrored in all 4 quadrants

    h=float(fft(fft(bzd,-1)/qq,1))
  aax0(0:nx-1,1:ny-1)= (h(0:nx-1,1:ny-1)-h(0:nx-1,0:ny-2))/dely
  if (ybc eq 1) then begin
    aax0(0:nx-1,     0)= (h(0:nx-1,     0)-h(0:nx-1,  ny-1))/dely   ; Periodic in y
    aax0(0:nx-1,    ny)=aax0(0:nx-1,0)
  endif
  
  aay0(1:nx-1,0:ny-1)=-(h(1:nx-1,0:ny-1)-h(0:nx-2,0:ny-1))/delx
  if (xbc eq 1) then begin
    aay0(     0,0:ny-1)=-(h(     0,0:ny-1)-h(  nx-1,0:ny-1))/delx   ; Periodic in x
    aay0(    nx,0:ny-1)=aay0(0,0:ny-1)
  endif 

  
 
    for j=0,ny-1 do begin
      aay0(*,j)=aay0(*,j)+imb(t)*(findgen(nx+1))*delx
    endfor

  
  
  writeu,1,aax0
  writeu,1,aay0
  
  bz_test=fltarr(nx,ny)
  bz_test=(aay0(1:nx,*)-aay0(0:nx-1,*))/delx - (aax0(*,1:ny)-aax0(*,0:ny-1))/dely
  wset,0
  tv,bytscl(bz_test,-100,100)
  print,'Frame ',strtrim(t,2),' processed'
  
endfor

close,1
print,'Photospheric boundary computed'






end
