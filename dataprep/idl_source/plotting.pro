pro plotting
@cleanup.blk
;This procedure computes (and plots) flux, flux imbalance, tilt angle and
;bipole separation as a function of time for the cleaned magnetograms (CLEANED)
;and the raw magnetograms (DATA). These plots are placed within the directory
;'plots'. The procedure also produces a series of stills comparing the cleaned
; and raw magnetograms for each frame. These are placed within the directory
;'frames'.

da=deltx^2 ; area of pixel in cm (to calculate flux)

dir2=dir+'plots/'

file_mkdir,dir2

;################# FLUX #########################################
print,''
print,'Producing plots:'


properties,data, pfluxdata,nfluxdata,imbdata,tiltd,sepd,cof_pd_x,cof_pd_y,cof_nd_x,cof_nd_y

properties,cleaned, pfluxclean,nfluxclean,imbclean,tiltc,sepc,cof_pc_x,cof_pc_y,cof_nc_x,cof_nc_y


fluxdata=(pfluxdata+nfluxdata)/2
fluxclean=(pfluxclean+nfluxclean)/2


set_plot,'ps' 
;Plot fluxes of cleaned magnetograms
device,filename=dir2+'cleaned_flux.eps',/encapsulated,/color,bits=8

plot,fluxclean/1.e20,xtitle='Frame',ytitle='Flux (10!U20!N Mx)',title='Cleaned Magnetogram',charsize=2,yrange=[0,1.3*max([nfluxclean,pfluxclean])/1.e20],/nodata
loadct,13,/silent
oplot,nfluxclean/1.e20,color=180
oplot,pfluxclean/1.e20,color=255
loadct,0,/silent
oplot,abs(imbclean)/1.e20,linestyle=2
al_legend,['Positive','Negative','|Imbalance|'],linestyle=[0,0,2],color=['red','green','black'],/right_legend
device,/close

print,dir2+'cleaned_flux.eps'




;Plot fluxes of raw magnetograms
device,filename=dir2+'raw_flux.eps',/encapsulated,/color,bits=8

plot,fluxdata/1.e20,xtitle='Frame',ytitle='Flux (10!U20!N Mx)',title='Raw Magnetogram',charsize=2,yrange=[0,max([nfluxdata,pfluxdata])*1.3/1.e20]
loadct,13,/silent
oplot,pfluxdata/1.e20,color=255
oplot,nfluxdata/1.e20,color=180
loadct,0,/silent
oplot,abs(imbdata)/1.e20,linestyle=2
al_legend,['Positive','Negative','|Imbalance|'],linestyle=[0,0,2],color=['red','green','black'],/right_legend
device,/close

print,dir2+'raw_flux.eps'






;Plot flux imbalance
device,filename=dir2+'flux_imbalance.eps',/encapsulated,/color,bits=8

!p.multi=[0,2,1]

range=max(abs([imbdata,imbclean]))/1.e20

plot,imbclean/1.e20,xtitle='Frame',ytitle='Flux Imbalance (10!U20!N Mx)',yrange=[-range,1.3*range],charsize=1
oplot,imbdata/1.e20,linestyle=2
al_legend,['Cleaned','Raw'],linestyle=[0,2],/right_legend


range=max(abs([imbdata/fluxdata,imbclean/fluxclean]))

plot,imbclean/fluxclean,xtitle='Frame',ytitle='Fractional Flux Imbalance',yrange=[-range,1.3*range],charsize=1
oplot,imbdata/fluxdata,linestyle=2
al_legend,['Cleaned','Raw'],linestyle=[0,2],/right_legend
device,/close
!p.multi=0

print,dir2+'flux_imbalance.eps'

set_plot,'x'






;plot tilt angle
set_plot,'ps'
device,filename=dir2+'tilt_angle.eps',/encapsulated
plot,tiltc/!dtor,xtitle='Frame',ytitle='Tilt Angle (Degrees)',yrange=[min([tiltc,tiltd]),1.3*max([tiltd,tiltc])]/!dtor,charsize=2
oplot,tiltd/!dtor,linestyle=2
al_legend,['Cleaned','Raw'],linestyle=[0,2],/right_legend
device,/close

print,dir2+'tilt_angle.eps'





;plot separation
device,filename=dir2+'separation.eps',/encapsulated
plot,sepc,xtitle='Frame',ytitle='Separation (Mm)',yrange=[0,1.3*max([sepd,sepc])],charsize=2
oplot,sepd,linestyle=2
al_legend,['Cleaned','Raw'],linestyle=[0,2],/right_legend
device,/close

print,dir2+'separation.eps'


set_plot,'x'









; PLOT AR PROPERTIES TO SCREEN
;print,'--------------------------------------------------------------'
print,''
print,'Window 5 displays the active region properties for the cleaned'
print,'and raw magnetograms.' 
print,''

!p.multi=[0,2,2]

window,5,xs=1200,ys=800
plot,fluxdata/1.e20,xtitle='Frame',ytitle='Flux (10!U20!N Mx)',title='Raw Magnetogram',charsize=2,yrange=[0,max([pfluxdata,nfluxdata])*1.3/1.e20],/nodata
loadct,13,/silent
oplot,pfluxdata/1.e20,color=255
oplot,nfluxdata/1.e20,color=180
loadct,0,/silent
oplot,abs(imbdata)/1.e20,linestyle=2
al_legend,['Positive','Negative','|Imbalance|'],linestyle=[0,0,2],color=['red','green','white'],/right_legend



plot,fluxclean/1.e20,xtitle='Frame',ytitle='Flux (10!U20!N Mx)',title='Cleaned Magnetogram',charsize=2,yrange=[0,max([fluxclean,nfluxclean,pfluxclean])/1.e20],/nodata
loadct,13,/silent
oplot,nfluxclean/1.e20,color=180
oplot,pfluxclean/1.e20,color=255
loadct,0,/silent
oplot,abs(imbclean)/1.e20,linestyle=2

al_legend,['Positive','Negative','|Imbalance|'],linestyle=[0,0,2],color=['red','green','white'],/right_legend


plot,tiltc/!dtor,xtitle='Frame',ytitle='Tilt Angle (Degrees)',yrange=[min([tiltc,tiltd]),1.3*max([tiltd,tiltc])]/!dtor,charsize=2
oplot,tiltd/!dtor,linestyle=2
al_legend,['cleaned','raw'],linestyle=[0,2],/right_legend

plot,sepc,xtitle='Frame',ytitle='Separation (Mm)',yrange=[0,1.3*max([sepd,sepc])],charsize=2
oplot,sepd,linestyle=2
al_legend,['cleaned','raw'],linestyle=[0,2],/right_legend


!p.multi=0






;produce .png images of the raw and cleaned frames

print,'--------------------------------------------------------------'
print,'Placing movie stills in ', dir+'frames'
window,0,xs=500,ys=250

file_mkdir,dir+'frames'


!p.background=255
!p.multi=[0,2,1]
for t=0,nt-1 do begin
num = strtrim(string(format='(I3.3)',t),2)
  pg_plotimage,bytscl(data(*,*,t),-100,100),color=0,/isotropic,title='Raw'
  loadct,13,/silent
  plots,cof_pd_x(t),cof_pd_y(t),psym=2,color=255
  plots,cof_nd_x(t),cof_nd_y(t),psym=2,color=150
  loadct,0,/silent
  pg_plotimage,bytscl(cleaned(*,*,t),-100,100),color=0,/isotropic,title='Cleaned'
  loadct,13,/silent
  plots,cof_pc_x(t),cof_pc_y(t),psym=2,color=255
  plots,cof_nc_x(t),cof_nc_y(t),psym=2,color=150
  loadct,0,/silent
  xyouts,5,5,'Frame '+num,color=0,/device,charsize=1
  write_png,dir+'frames/frame'+num+'.png',tvrd(/true)
endfor
!p.background=0
!p.multi=0


print,''



end


