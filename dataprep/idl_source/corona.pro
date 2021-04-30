pro corona
;previews _00000p file produced by corona.f90
opt=1l
alpha=0.d0
dir='--------------------'

openr,10,'params.dat'
readf,10,nmax,nt,alpha,dir
close,10


nsize=2*2^nmax



aax=fltarr(nsize,nsize+1,nsize+1)
aay=fltarr(nsize+1,nsize,nsize+1)
aaz=fltarr(nsize+1,nsize+1,nsize)

filename='run1_00000p'
filename=strtrim(dir,2)+'hexa_files/run1_00000p'
print,filename
openr,3,filename,/f77_unformatted

readu,3,opt
readu,3,aax
readu,3,aay
readu,3,aaz
close,3

dx=6./nsize

window,1,xs=nsize,ys=nsize
;restore,'mag_data.sav'
;restore,'prepot.sav'
bz = dblarr(nsize,nsize,nsize)
for k=0,10 do begin
for i=0,nsize-1 do begin
  for j=0,nsize-1 do begin
    bz(i,j,k) = (aay(i+1,j,k)-aay(i,j,k))/dx - (aax(i,j+1,k) - aax(i,j,k))/dx
  endfor
endfor
tv,bytscl(bz(*,*,k),-100,100)
wait,0.25
endfor


stop
end