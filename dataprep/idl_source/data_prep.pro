pro data_prep
@cleanup.blk
;This routine takes an IDL save (data.sav) file containing an array called data
;which is fltarr(nx,ny,nt) containing de-rotated magnetograms and prepares them
;to be used as the lower boundary condition in Hexa. 
;
;Firstly, the magnetograms are cleaned up. Various options are available
;
;Secondly various plots are produced comparing the properties of the raw
;magnetograms with the cleaned magnetograms. These are placed in the 'plots/'
;sub-directory. Additionally .PNG images are produced comparing each raw and cleaned
;frame. These are placed in the 'frames/' directory.
;
;Thirdly, the cleaned magnetograms are prepared for use in Hexa. This involves
;placing them into an appropriately-sized computational domain (128^2, 256^2 or 512^2)
;then flux balancing is applied. Various options are available. Also, the user is
;given the opportunity to set the start and end frames for the simulation.
;
;Lastly the file 'mag_data' is written out for use in the fortran codes
;'evolve.f90' and 'corona.f90' which produce the '_evol' file and a non-linear
;force-free initial condition for use with Hexa. These fortran codes are run by
;this IDL routine, and the outputs are placed in the 'hexa_files' subdirectory.
;
;Finally, plots of the properties of the lower boundary condition are produced
;and are placed in the plots subdirectory

set_plot,'x'
device,decomposed=0
loadct,0,/silent
close,/all



print,''
print,'##############################################################'
print,'#     Welcome to the clean-up and data preparation code!     #'
print,'##############################################################'
print,''
@required_info
print,'The working directory is: ',dir
print,'reading in '+dir+'data.sav...'
restore,dir+'data.sav' ;IDL save file containing the variable 'DATA' which is a fltarr(nx,ny,nz)


datainfo=size(data)

nx=datainfo(1) ; get nx
ny=datainfo(2) ; get ny
nt=datainfo(3) ; get nt

print,'Magnetogram data contains (nx,ny,nz) = (',strtrim(nx,2),',',strtrim(ny,2),',',strtrim(nt,2),')'

cleaned=fltarr(nx,ny,nt) ;variable to store cleaned magnetograms
;cleaned=data




;preview the data - if the user chooses to
response='a'
print,''
jump1:print,'Would you like to preview the data? (y/n)'
read,response
if response eq 'y' then begin
   preview ;call preview procedure (contained in preview.pro)
endif else if response ne 'n' then begin
   print,'Invalid response! Please try again...'
   goto,jump1
endif



jump2:print,''



print,'##################### Cleaning Procedure #####################'

cleaning ; call cleaning procedure (contained within cleaning.pro)









print,'##################### Plotting Procedure #####################'

plotting ; call plotting procedure (contained within plotting.pro)
answer='a'
print,'--------------------------------------------------------------'
Jump3: print,'Would you like to re-run the cleaning? (y/n)'
read,answer
if answer eq 'y' then begin
   wdelete,0
   goto,jump2 
endif else if answer ne 'n' then begin
   print,'Invalid response! Please try again...'
   goto,jump3
endif
wdelete,0

print,''




print,'###################### Data preparation ######################'

;stop

preparation ; call preparation procedure (contained within preparation.pro)

print,'--------------------------------------------------------------'
print,''


wdelete,5

;reset nx,ny to hexa grid values, and define nmax (needed by fortran codes) 
if size eq 128 then begin
    nmax=6
    nx=128
    ny=128
endif else if size eq 256 then begin
    nmax=7
    nx=256
    ny=256
endif else begin
    nmax=8
    nx=512
    ny=512
endelse

if open ne 1 and periodic ne 1 then begin

;select alpha
print,'An initial linear force-free coronal field must now be generated.'
alpha1:print,'What value of alpha would you like to use? (range=[-3.3,3.3])'
read,alpha
if alpha gt 3.3 or alpha lt -3.3 then begin
  print, 'Invalid alpha chosen! Please try again...
  goto,alpha1
endif

endif else alpha = 0.0

file_mkdir,dir+'hexa_files'

;write out parameter file read in by fortran codes 
openw,10,'params.dat'
printf,10,nmax,nt,alpha,'  ',dir,'                   '
printf,10,periodic,open
close,10

;write out magnetogram data to fortran readable file
openw,10,dir+'mag_data.dat',/f77
writeu,10,double(cor_mag)
close,10


if open eq 1 or periodic eq 1 then begin
  ; run idl fourier_method.pro
  print,'--------------------------------------------------------------'
  print,'Now running fourier_method.pro to generate the Hexa files:'
  print,''
  len=strlen(dir)
  len=len-1
  file=strmid(dir,0,len)
  
  fourier_method,cor_mag,dir,file,periodic,open
  print,'--------------------------------------------------------------'
endif else begin

    ;run fortran codes
    print,'--------------------------------------------------------------'
    print,'Now running the fortran codes to generate the Hexa files:'
    print,''
    spawn,'./evolve'
    print,'--------------------------------------------------------------'
    print,''
    spawn,'./corona'
    print,''
    print,'--------------------------------------------------------------'

endelse

;write param1
openw,10,dir+'hexa_files/param1'
printf,10,2l
printf,10,size, size, size
printf,10,l0
printf,10,dtime*60.
close,10

;write _setup file
len=strlen(dir)
openw,10,dir+'hexa_files/'+strmid(dir,0,len-1)+'_setup'
printf,10,'vsetup=3'
;printf,10,'nrelax=0'
printf,10,'nmajor=',strtrim(nt,2)

vsetup=3
nmajor=120
printf,10,'nstrt=0'
printf,10,'nend=',strtrim(nt-1,2)
printf,10,'etaia=0.00'
printf,10,'eta4a=0.000'
printf,10,'periodic=',strtrim(periodic,2)
printf,10,'open=',strtrim(open,2)


;printf,10,'nminor=1'
;printf,10,'nrep=500'
;printf,10,'nsplit=0.0'
;printf,10,'nbpole=0.0'
;printf,10,'orient=180'
;printf,10,'imode=2'
;printf,10,'nsinit=2'
;printf,10,'seed=0'
;printf,10,'rad=0.6'
;printf,10,'nfil=0.0'
;printf,10,'bxcon=1.0'
close,10


Print,'Written param1 and ',+strmid(dir,0,len-1)+'_setup into ',dir+'hexa_files/'


print,''

;plot properties of AR to screen
deltx=l0*6./nx

properties,cor_mag, pflux,nflux,imb,tilt,sep,cof_pd_x,cof_pd_y,cof_nd_x,cof_nd_y

window,5,xs=1600,ys=500
!p.multi=[0,3,1]
plot,pflux/1.e20,xtitle='Frame',ytitle='Flux (10!U20!N Mx)',charsize=2,yrange=[0,max(pflux/1.e20)]
plot,tilt/!dtor,xtitle='Frame',ytitle='Tilt angle (degrees)',charsize=2
plot,sep,yrange=[0,max(sep)],xtitle='Frame',ytitle='Separation (Mm)',charsize=2
!p.multi=0



print,'##########################  Summary ##########################'

;and plot to file
set_plot,'ps'
print,'producing plots:'

device,filename=dir+'plots/hexa_flux.eps',/encapsulated

plot,pflux/1.e20,xtitle='Frame',ytitle='Flux (10!U20!N Mx)',charsize=2,yrange=[0,max(pflux/1.e20)]
print,dir+'plots/hexa_flux.eps'

device,/close
device,filename=dir+'plots/hexa_tilt.eps',/encapsulated
plot,tilt/!dtor,xtitle='Frame',ytitle='Tilt angle (degrees)',charsize=2
print,dir+'plots/hexa_tilt.eps'

device,/close
device,filename=dir+'plots/hexa_sep.eps',/encapsulated
plot,sep,yrange=[0,max(sep)],xtitle='Frame',ytitle='Separation (Mm)',charsize=2
print,dir+'plots/hexa_sep.eps'
device,/close
set_plot,'x'

file_copy,'params.dat',dir+'hexa_files/',/overwrite

print,''
print,'Hexa length unit (cm): ',strtrim(l0,2)
print,'Number of frames to be used in Hexa: ', strtrim(nt,2)
print,'Hexa resolution: ',strtrim(size,2),'^3'
print,'Hexa files located in '+dir+'hexa_files/'
print,''
if (open eq 1) or (periodic eq 1) then begin
   print,' NB: Please use the codes in potenFourier/ to produce the potential fields'
endif else begin
   print,' NB: Please use the codes in potenMPI/ to produce the potential fields'
endelse
print,'Done!'

;stop
end










