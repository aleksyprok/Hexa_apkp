pro properties, array,pflux,nflux,imb,tilt,sep,pposx,pposy,nposx,nposy
;calculate the flux, tile angles and separation of a time series of magnetograms
@cleanup.blk
dum=size(data)


pflux=fltarr(nt)
nflux=fltarr(nt)
imb=fltarr(nt)
tilt=fltarr(nt)
sep=fltarr(nt)
pposx=fltarr(nt)
pposy=fltarr(nt)
nposx=fltarr(nt)
nposy=fltarr(nt)


da=deltx^2
;flux imbalance
for t=0,nt-1 do begin
   imb(t)=total(array(*,*,t))*da
   pflux(t)=total(abs(array(*,*,t)))*da
endfor

;positive and negative flux
pflux=pflux/2.
nflux=pflux

pflux=pflux+imb/2
nflux=nflux-imb/2


for t=0,nt-1 do begin
;print,'t=',t
  dta=array(*,*,t)
  
  
  dtap=dta
  dtap(where (dta lt 0)) = 0
  dtan=dta
  dtan(where (dta gt 0)) = 0

 
;calculate centre of flux of positive polarity in raw magnetograms
  cof,dtap,nx,ny,x,y 
  pposx(t)=x
  pposy(t)=y

;calculate centre of flux of negative polarity in raw magnetograms
  cof,dtan,nx,ny,x,y
  nposx(t)=x
  nposy(t)=y

endfor

;calculate tilt angle
tilt=atan( (nposy - pposy ) , (  nposx - pposx   ) ) ; raw 

;separation
sep= sqrt( (nposy - pposy )^2 + (nposx - pposx)^2)*deltx/1.e8 ;raw

end