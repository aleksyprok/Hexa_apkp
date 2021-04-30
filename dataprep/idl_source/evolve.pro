pro evolve
;previews _evolve file produced by 'evolve.f90'
nmax=0
nt=0
alpha=0.d0
dir='--------------------'

openr,10,'params.dat'
readf,10,nmax,nt,alpha,dir
close,10


nsize=2*2^nmax

print,nsize,nt

aax=fltarr(nsize,nsize+1)
aay=fltarr(nsize+1,nsize)


window,1,xs=nsize,ys=nsize


fname=strtrim(dir,2)+'hexa_files/run1_evolve'
print,fname
opt=0l
openr,10,fname,/f77
readu,10,opt

for t=0,nt-1 do begin

readu,10,aax
readu,10,aay




dx=6./float(nsize)

bz=fltarr(nsize,nsize)

;help,aax,aay,bz
bz=(aay(1:nsize,*)-aay(0:nsize-1,*))/dx - (aax(*,1:nsize)-aax(*,0:nsize-1))/dx
print,t,max(bz)

tv,bytscl(bz,-100,100)

wait,0.05

endfor
close,10

end