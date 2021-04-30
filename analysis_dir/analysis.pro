pro analysis

dir ='../hexaf90/run1' ;directory containing hexa files

cadence=1.
l0=1.

root='param1'

openr,10,dir+'/'+root
dum=1
readf,10,dum
readf,10,dum,dum,dum
readf,10,l0
readf,10,cadence
close,10

l0=l0/1.e8 ; convert to Mm
cadence=cadence/60 ;convert to minutes

;Plot Helicity

root = 'helicity'
filename=dir+'/'+root

openr,10,filename,/f77_unformatted
n = long(0)
readu,10,n
helicity = dblarr(n)
readu,10,helicity
close,10





tt=findgen(n)/(n-1)
tt=tt*cadence*n/60./24.



set_plot,'ps'
device,filename=dir+'/Helicity.eps',/encapsulated,xsize=7,ysize=5,/inches

plot,tt,helicity,xtitle='Time (Days)',ytitle='Helicity (Mx!U2!n)',charsize=1.5
device,/close
set_plot,'x'

fe=dblarr(n)
te=dblarr(n)

dum=0
nt=0
openr,10,dir+'/params.dat'
readf,10,dum,nt
close,10

start = 0
stop =nt-1




for i=start,stop do begin

  ;read in energy files
  root = 'fenergy'
  ending = 'p'
  num = strtrim(string(format='(I5.5)',i),2)
  filename=dir+'/'+root+'_'+num+ending

  nx = long(0)
  ny = long(0)
  nz = long(0)

  openr,10,filename,/f77_unformatted
  readu,10,nx,ny,nz

  fenergy = dblarr(nx+1,ny+1,nz+1)
  b2 = dblarr(nx+1,ny+1,nz+1)
  cc=dblarr(nx+1,ny+1,nz+1)

  readu,10,fenergy
  readu,10,b2
  readu,10,cc


  close,10

  ; read in 'base' which contains the magnetograms
  if i eq start then begin
  root = 'base'
  fname = dir+'/'+root
  openr,10,fname,/f77_unformatted
  n=long(0)
  readu,10,n
  base=dblarr(nx+1,ny+1,n)
  readu,10,base
  close,10
  endif


  ;make map of the current density^2

  xx=findgen(nx+1)/float(nx)*6.*l0
  yy=findgen(ny+1)/float(ny)*6.*l0

  pic = total(cc(*,*,*),3)
  ;pic = pic+1

  ;pic = alog(pic) ;log image to bring out detail
  print,max(pic),min(pic)
  up=max(pic)/2
  down=0.


  set_plot,'ps'
  root = 'current'
  ending = '.eps'
  num = strtrim(string(format='(I5.5)',i),2)
  filename=dir+'/'+root+'_'+num+ending

  device,/encapsulated,/color,bits=8,xsize=5,ysize=5,filename=filename,/inches

  loadct,3,/silent

  pg_plotimage,255-bytscl(pic,down,up),xx,yy,/isotropic,background=255,thick=4,xthick=4,ythick=4,xtitle='x (Mm)',ytitle='y (Mm)'

  loadct,13,/silent
  contour,base(*,*,i),xx,yy,levels=[-411,-205,-102,-51,-25.6],/overplot,/isotropic,color=70,c_thick=2
  contour,base(*,*,i),xx,yy,levels=[25.6,51,102,205,411],/overplot,/isotropic,color=150,c_thick=2
loadct,0,/silent




  device,/close

 ;make map of free magnetic energy storage

  root = 'fenergy'
  ending = '.eps'
  num = strtrim(string(format='(I5.5)',i),2)
  filename=dir+'/'+root+'_'+num+ending

  device,filename=filename,/color,bits=8

  pic=total(fenergy,3)

  up=mean(abs(pic))+stddev(pic);stddev(pic);max(pic)
  down=-mean(abs(pic))-stddev(pic);-stddev(pic);min(pic)

  pg_plotimage,bytscl(pic,down,up),xx,yy,/isotropic,xtitle='x (Mm)',ytitle='y (Mm)'
  loadct,13,/silent
  contour,base(*,*,i),xx,yy,levels=[-411,-205,-102,-51,-25.6],/overplot,/isotropic,color=70,c_thick=2
  contour,base(*,*,i),xx,yy,levels=[25.6,51,102,205,411],/overplot,/isotropic,color=150,c_thick=2
loadct,0,/silent

  device,/close
  set_plot,'x'

  fe(i) = total(fenergy)
  print,i

  te(i) = total(b2)


endfor

loadct,0,/silent

set_plot,'ps'
device,filename=dir+'/Free_Energy.eps',/encapsulated,xsize=7,ysize=5,/inches
plot,tt,fe,xtitle='Time (Days)',ytitle='Free Magnetic Energy (ergs)',charsize=1.5
device,/close
set_plot,'x'




stop
end
