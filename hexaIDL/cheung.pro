pro cheung,nlines=nlines ;nlines is by default 10,000. 

;Produces emission proxy images using the method described by Cheung & DeRosa 2012.

; the array 'emission' is produced, which is a 3d array containing the emission as a function of position throughout the volume. By default this array is integrated along the z direction to obtain the emission from above, but it can in principle be integrated in any direction to obtain the emission from any line of sight.

@hexa.blk
@hexa.trc
common chng, emission,upres

upres=1 ;this parameter sets the resolution of the image produced relative to the simulation.
; if upres=1 then the image will have the same resolution as the simulation (e.g. if the
; simulation was 256x256x256 then the image will be 256x256). If upres=2 then if the simulation
; had a resolution of 256x256x256 then the imsge will have a resolution of 512x512.
;essentially use this option if you want a prettier looking image for publication

neach=100 ;number of field lines calculated between renderings



emission=fltarr((nx+1)*upres,(ny+1)*upres,(nz+1)*upres)
window,1,xs=600,ys=600


if keyword_set(nlines) then begin
   nits=nlines/neach
   if nits lt 1 then begin
      print,'error, too few fieldlines chosen!'
      stop
   endif
endif else nits=100



seed=3141592l
counter=0l



for dum=0,nits-1 do begin ;loop over iterations

    xs=randomu(seed,neach)*6 ;randomly pick points in x
    ys=randomu(seed,neach)*6 ;randomly pick points in y

    for i=0,neach-1 do begin ;loop over neach
      x=xs(i) ;x startpoint for field line
      y=ys(i) ;y startpoint for field line
      z=delx  ;z startpoint for fieldline (picked to be 1 gridcell up from the photosphere)
      
      counter=counter+1l
      dem,x,y,z,dq,upres ;get the emission from one field line (dq)
      emission=emission+dq ;add the emission from one field line to the total emission
    
    endfor
    
    ;preview image from this iteration after neach fieldlines have been added to the emission 
    
    pic=(total(emission,3)+1) ;integrate in z direction
    pic=alog(pic) ;log image

    up=max(pic)/1.
    down=0.


    device,decomposed=0
    loadct,3,/silent

    xx=findgen((nx+1)*upres)*delx/upres
    yy=findgen((ny+1)*upres)*dely/upres

    pg_plotimage,bytscl((pic),down,up),xx,yy
    loadct,13,/silent
    xx=findgen(nx+1)*delx
    yy=findgen(ny+1)*dely
    contour,bz1(*,*,0),xx,yy,levels=[25,50,100,200],color=180,/overplot
    contour,bz1(*,*,0),xx,yy,levels=[-200,-100,-50,-25],color=100,/overplot

    loadct,0,/silent
    device,decomposed=1
    print,counter ;print total number of field lines used so far

endfor


;display final image


bb=sqrt(bx1^2+by1^2+bz1^2)
bb=congrid(bb,(nx+1)*upres,(ny+1)*upres,(nz+1)*upres)

pic=(total(emission*bb,3)+1)
pic=total(emission,3)+1
;pic=alog(pic)

up=max(pic)/4.
down=0.

;up=max(pic)
;down=max(pic)/2

!p.background=255
window,1,xs=600,ys=600
device,decomposed=0
loadct,3,/silent

xx=findgen((nx+1)*upres)*delx/upres
yy=findgen((ny+1)*upres)*dely/upres

pg_plotimage,255-bytscl((pic),down,up),xx,yy,color=0
loadct,13,/silent
xx=findgen(nx+1)*delx
yy=findgen(ny+1)*dely
contour,bz1(*,*,0),xx,yy,levels=[25,50,100,200],color=180,/overplot
contour,bz1(*,*,0),xx,yy,levels=[-200,-100,-50,-25],color=100,/overplot

loadct,0,/silent
device,decomposed=1


end






pro dem,x,y,z,dq,upres
@hexa.trc
@hexa.blk
;calculates the emission along a field line (dq) plotted from x,y,z


trace,x,y,z ;trace field line from (x,y,z)

n=n_elements(xtrc)

if ztrc(n-1) lt 1 then begin ;if the end point of the field line is the photosphere (closed field line)

    cx=fltarr(n) ;j_x along field line
    cy=fltarr(n) ;j_y along field line
    cz=fltarr(n) ;j_z along field line


    for i=0,n-1 do begin ;move along field line
      tr_magn,xtrc(i),ytrc(i),ztrc(i),dum,cx2,cy2,cz2,par=1 ;interpolate current on field line
      cx(i)=cx2
      cy(i)=cy2
      cz(i)=cz2
    endfor

    s=total(ds2,/cumulative) ;length of field line

    c2=cx^2+cy^2+cz^2 ;j^2


    meanj=total(c2*ds2)/total(ds2) ;mean current along field line

    dq=fltarr((nx+1)*upres,(ny+1)*upres,(nz+1)*upres)*0.0

    for i=0,n-1 do begin ;loop along field line
	xc=floor(xtrc(i)/delx*upres) ;x-index of this point along the field line in dq array
	yc=floor(ytrc(i)/dely*upres) ;y-index...
	zc=floor(ztrc(i)/delz*upres) ;z-index...

	dq(xc,yc,zc)=meanj ;set this location in dq to be the mean j emission
    endfor

endif else dq=0. ;open field lines produce no emission (as per Cheung & DeRosa)

end