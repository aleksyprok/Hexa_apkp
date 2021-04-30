pro preview
@cleanup.blk
;previews magnetogram data. Graphs the flux, and shows a movie of the magnetograms


print,'Previewing magnetogram data...'
print,''
print,"Window 0 displays the active region's flux."
print,''
print,'Window 1 displays a movie of the magnetograms.'
print,''

;calculate flux
flux=fltarr(nt)
imb=fltarr(nt)

for t=0,nt-1 do begin
    flux(t) = total(abs(data(*,*,t)))*deltx^2/2.
    imb(t) = total(data(*,*,t))*deltx^2
endfor

flux=flux/1.e21
imb=imb/1.e21

;plot flux to screen
window,0,xs=800,ys=500
plot,flux,xtitle='Frame',ytitle='Flux (10!U21!N Mx)',charsize=2,yrange=[0,max([flux+imb/2,flux-imb/2])],/nodata
loadct,13,/silent
oplot,flux+imb/2,color=255
oplot,flux-imb/2,color=180
loadct,0,/silent
oplot,abs(imb),linestyle=2
al_legend,['Positive','Negative','|Imbalance|'],linestyle=[0,0,2],color=['red','green','white'],/right_legend



;movie of magnetograms
!p.background=255
window,1,xs=500,ys=500

xx=findgen(nx)*deltx/1.E8
yy=findgen(ny)*deltx/1.e8

jump1:for t=0,nt-1 do begin
   pg_plotimage,bytscl(data(*,*,t),-100,100),xx,yy,xtitle='x (Mm)',ytitle='y (Mm)',title='Frame '+strtrim(t+1,2)+' of '+strtrim(nt,2),color=0,/isotropic
   wait,0.05
endfor
response='a'

jump2:print,'Would you like to replay the movie? (y/n)'
read,response
if response eq 'y' then goto,jump1 else if response ne 'n' then begin
  print,'Invalid response! Please try again...'
  goto,jump2
endif
!p.background=0

jump3: print,'Would you like to continue? (y/n)'
read,response
if response eq 'n' then stop else if response ne 'y' then begin
    print,'Invalid response! Please try again...'
    goto,jump3
endif

wdelete,0,1


end

