;---------------------------------------------------------------
;  Evolution on hexagonal grid (HEXA) Main Menu:
;---------------------------------------------------------------
pro HEXA
common wid,wid_main
@hexa.blk
;
;  Create Main Menu widget:
;
base=widget_base(title='HEXA '+version,/column,xoffset=100,yoffset=50)
wid_main=lonarr(4)
wid_main(0)=base
dummy={MENU, flags: 0, name: '', proc: ''}
desc=[{MENU,1,'         Main Menu         '      ,'MENU'            },  $
      {MENU,0,'Set HEXA Home Directory'          ,'PROG_HOME'       },  $
      {MENU,0,'Setup Directory and Grid'         ,'PROG_GRID'       },  $
      {MENU,0,'Restore 3D Model from File'       ,'PROG_RESTORE'    },  $
      {MENU,0,'Restore Current File'             ,'PROG_RESTORE_CUR'},  $
      {MENU,0,'Recompute Magnetic Field'         ,'PROG_FIELD'      },  $
      {MENU,0,'Show 1D Display'                  ,'PROG_PLOT1D'     },  $
      {MENU,0,'Show 2D Display'                  ,'PROG_PLOT2D'     },  $
      {MENU,0,'Show 3D Display'                  ,'PROG_PLOT3D'     },  $
      {MENU,0,'Show Surface Evolution'           ,'PROG_EVOLVE'     },  $
      {MENU,0,'Save Field Lines'                 ,'PROG_SAVE_LINES' },  $
      {MENU,0,'HEXA Colour Table'                 ,'hexa_colour'      },  $
      {MENU,0,'About HEXA'                       ,'HEXA_ABOUT'      },  $
      {MENU,2,'Quit'                             ,'QUIT'            }]
menu=cw_pdmenu(base,desc,uvalue=desc.proc)
;
txt='Home:                         '       ; string with appropriate length
if n_elements(home) gt 0 then strput,txt,home,6
wid_main(1)=widget_label(base,value=txt,/align_left)
;
txt='Dir:                          '       ; string with appropriate length
if n_elements(dir) gt 0 then strput,txt,dir,5
wid_main(2)=widget_label(base,value=txt,/align_left)
;
txt='File:                         '       ; string with appropriate length
if n_elements(filename) eq 0 then filename=''
if filename ne '' then begin
  if n_elements(nt) eq 0 then nt=0L
  file=string(filename,nt,format='(a,"_",i5.5)')
  strput,txt,file,6
endif
wid_main(3)=widget_label(base,value=txt,/align_left)
widget_control,base,/realize
xmanager,'hexa',base,/no_block
end

;-------------------------------------------------------------------------
;  Process Main Menu Events:
;-------------------------------------------------------------------------
pro HEXA_EVENT,ev
common wid,wid_main
@hexa.blk
;
;  Get UVALUE of event. For pull-down menus UVALUE is an element
;  of a string array:
;
widget_control,ev.id,get_uvalue=uvalue
if strmid(uvalue(0),0,4) eq 'MENU' then uvalue=uvalue(ev.value)
;
;  Special cases for which no parameters need to be passed to
;  any procedure:
;
if uvalue eq 'HEXA_ABOUT' then begin
  title='Evolution on hexagonal grid (HEXA)'
  text=strarr(17)
  text( 0)='Evolution on hexagonal grid (HEXA)'
  text( 1)='Version 1.1 -- January 2008'
  text( 2)='Developed by:'
  text( 3)='  Aad van Ballegooijen'
  text( 4)='  Smithsonian Astrophysical Observatory'
  text( 5)='  60 Garden Street, MS 15'
  text( 6)='  Cambridge, MA 02138'
  text( 7)='  Email: vanballe@cfa.harvard.edu'
  text( 8)=''
  text( 9)='Version 2.0 -- May 2015'
  text(10)='Revised by:'
  text(11)='  Gordon Gibb'
  text(12)='  University of St Andrews'
  text(13)='  Mathematical Institute'
  text(14)='  St Andrews'
  text(15)='  KY16 9SS'
  text(16)='  Email: gpsg@st-andrews.ac.uk'



  POP_UP, title, text
  return
endif
if uvalue eq 'QUIT' then begin
  widget_control,ev.top,/destroy
  return
endif
;
;  All other programs:
;
if strmid(uvalue,0,5) eq 'PROG_' then begin
  widget_control,ev.top,/hourglass
  call_procedure,uvalue
  if uvalue eq 'PROG_HOME' then begin
    txt='Home:                         '  ; string with appropriate length
    if n_elements(home) gt 0 then strput,txt,home,6
    widget_control,wid_main(1),set_value=txt
  endif
  if uvalue eq 'PROG_GRID' then begin
    txt='Dir:                          '  ; string with appropriate length
    if n_elements(dir) gt 0 then strput,txt,dir,5
    widget_control,wid_main(2),set_value=txt
  endif
  txt='File:                         '  ; string with appropriate length
  if n_elements(filename) gt 0 and n_elements(nt) gt 0 then begin
    if filename ne '' then file=string(filename,nt,format='(a,"_",i5.5)') $
                      else file=''
    strput,txt,file,6
  endif
  widget_control,wid_main(3),set_value=txt
endif
;
end

;-------------------------------------------------------------------------
;  Display message in pop-up window:
;-------------------------------------------------------------------------
pro POP_UP,title,text
base=widget_base(title=title,/column,xoffset=300,yoffset=270)
if title eq 'Error' then begin
  wlis=widget_list(base,value=text,ysize=n_elements(text),  $
                        resource_name='red')
endif else begin
  wlis=widget_list(base,value=text,ysize=n_elements(text),  $
                        resource_name='white')
endelse
done=widget_button(base,value='OK',resource_name='magenta')
widget_control,base,set_uvalue=done,/realize
xmanager,'POP_UP',base,/no_block
end

;-------------------------------------------------------------------------
;  Event handler for pop-up window:
;-------------------------------------------------------------------------
pro POP_UP_EVENT,ev
widget_control,ev.top,get_uvalue=done
if ev.id eq done then widget_control,ev.top,/destroy
end

;------------------------------------------------------------------
;  Initialize colors:
;------------------------------------------------------------------
pro hexa_colour
@hexa.col
;
;  For PostScript device:
;
if !d.name eq 'PS' then begin
  r=bytarr(256)
  g=bytarr(256)
  b=bytarr(256)
  nfbl=10          ; number of colors for blue scale (dips)
;
  ncol=256-nfbl-6
  fgrey=findgen(ncol)/(ncol-1)
  r(0:ncol-1)=round(255.*fgrey)    ; grey scale for images
  g(0:ncol-1)=r(0:ncol-1)
  b(0:ncol-1)=r(0:ncol-1)
  fblue=findgen(nfbl)/(nfbl-1)
  r(ncol:255)=[intarr(nfbl)     ,255,  0,  0,  0,255,255]
  g(ncol:255)=[round(255.*fblue),  0,255,  0,255,  0,255]
  b(ncol:255)=[intarr(nfbl)+255 ,  0,  0,255,255,255,  0]
  tvlct,r,g,b
  nsum=ncol+nfbl
  black  =0
  white  =nsum-1
  red    =nsum+0
  green  =nsum+1
  blue   =nsum+2
  ltblue =nsum+3
  magenta=nsum+4
  yellow =nsum+5
endif
;
;  For X device (always use TrueColor Visual):
;
if !d.name eq 'X' then begin
;
  device,true_color=24,decomposed=1   ; do before any window is opened
  loadct,0,/silent                    ; use ONLY this color table
  nfbl=0
  ncol=256
  black  =  0+256L*(  0+256L*  0)     ; special colors
  white  =255+256L*(255+256L*255)
  red    =255+256L*(  0+256L*  0)
  green  =  0+256L*(255+256L*  0)
  blue   =  0+256L*(  0+256L*255)
  ltblue =  0+256L*(255+256L*255)
  magenta=255+256L*(  0+256L*255)
  yellow =255+256L*(255+256L*  0)
;
endif
;
;  Color names (not needed):
;
;col_names=['Black','White','Red','Green','Blue','LtBlue','Magenta','Yellow']
;col_index=[ black , white , red , green , blue , ltblue , magenta , yellow ]
;
;  Field-line colors:
;
colors=[blue,ltblue,magenta]
end

;------------------------------------------------------------------
;  Apply grey scale:
;------------------------------------------------------------------
function grey_scale,a,amin,amax,gamma=gamma,neg=neg
@hexa.col
;
;  Usage:  tv,grey_scale(a,amin,amax),....
;
if n_params() eq 1 then amax=max(a,min=amin)
if n_elements(neg) eq 0 then neg=0
if neg lt 2 then begin
  frac=(float(a-amin)/(amax-amin))>0.<1.
  if n_elements(gamma) gt 0 then frac=frac^gamma
  if neg eq 1 then frac=1.0-frac
endif else begin
  amin1=float(amin)>0.001
  amax1=float(amax)>(amin1+0.001)
  frac=(alog((float(a)>0.001)/amin1)/alog(amax1/amin1))>0.<1.
endelse
return,byte((ncol-1)*frac)
end

;------------------------------------------------------------------
;  Apply color scale for SQL maps:
;------------------------------------------------------------------
function color_scale,a
@hexa.col
;
;  Usage:  tv,color_scale(a),true=1,....
;
siz=size(a)
nx=siz(1)
ny=siz(2)
img=bytarr(3,nx,ny)                                   ; black
;
amax=max(a)
;
for j=0,ny-1 do begin
  for i=0,nx-1 do begin
    if a(i,j) eq -3.0 then img(*,i,j)=[  0,150,  0]   ; green  (top)
    if a(i,j) eq -2.0 then img(*,i,j)=[100,100,  0]   ; orange (side)
    if a(i,j) eq -1.0 then img(*,i,j)=[255,255,255]   ; white  (short)
    if a(i,j) gt 0.0 then begin
      frac=(a(i,j)/amax)^0.1
      img(2,i,j)=byte(fix(ncol*1.5*frac     )>0<(ncol-1))
      img(1,i,j)=byte(fix(ncol*(2.*frac-1.0))>0<(ncol-1))
      img(0,i,j)=byte(fix(ncol*(4.*frac-3.0))>0<(ncol-1))
    endif
  endfor
endfor
;
return,img
end

;----------------------------------------------------------------
;  Set home directory:
;----------------------------------------------------------------
pro prog_home
@hexa.blk
;
;  Set home directory:
;
if n_elements(home) gt 0 then default=home else cd,current=default
print,'Default directory:  '+default
home1=''
read,prompt='Enter directory:  ',home1
if home1 eq '' then home1=default
if strmid(home1,strlen(home1)-1,1) ne '/' then home1=home1+'/'
cd,home1
home=home1   ; change HOME only if CD was succesful
end

;-----------------------------------------------------------------
;  Setup hexagonal grid:
;-----------------------------------------------------------------
pro prog_grid,current=current
@hexa.blk
;
;  Select subdirectory and find grid parameter file (param1):
;
if n_elements(current) eq 0 then current=0
if current then begin
  dir1=dir
endif else begin
  if n_elements(dir) gt 0 then default=dir else default='Filament2'
  print,'Default directory:  '+default
  dir1=''
  read,prompt='Select directory name:  ',dir1
  if dir1 eq '' then dir1=default
  if strmid(dir1,strlen(dir1)-1,1) ne '/' then dir1=dir1+'/'
endelse
file1=dir1+'param1'
openr,1,file1,error=err
if err ne 0 then begin
  print,!err_string
  print,'Error: Cannot find '+file1
  close,1
  print,'-----'
  return
endif
dir=dir1
print,'Directory:     '+dir
;
;  Read grid parameter(s):
;
print,'Read grid parameter(s): ',file1
num_hex_cells=0
readf,1,num_hex_cells
nx=0L  &  ny=0L  &  nz=0L
readf,1,nx,ny,nz
close,1
print,string(num_hex_cells,format='("num_hex_cells=",i3)')
print,string(nx,ny,nz,format='("nx=",i3,", ny=",i3,", nz=",i3)')
;
xoff=0.5   ; offsets necessary for periodic boundary conditions
yoff=0.5
;
;  Position arrays and spacings in 3D grid:
;
imax=num_hex_cells
jmax=round(imax*sqrt(3.0))
xmin=0.0
xmax=3.0*imax-1.0+2.*xoff
;x2min=-xmax
;x2max=2.*xmax
x2min=xmin
x2max=xmax
ymin=0.0
ymax=xmax
dy=(ymax-2.*yoff)/float(2*jmax)
delx=xmax/nx
dely=ymax/ny
delz=delx
zmin=0.0
zmax=nz*delz
;
yc1=dy*findgen(2*jmax+1)+yoff
xc1=fltarr(2*jmax+1)
xc1(2*indgen(jmax+1))=0.5
x0=xoff
yc=[yc1,yc1]
xc=[x0+xc1,x0+2.0-xc1]
for i=1,imax-1 do begin
  x0=xoff+3.0*i
  yc=[yc,yc1,yc1]
  xc=[xc,x0+xc1,x0+2.0-xc1]
endfor
nc=imax*(4*jmax+2)
print,string(nc,format='("Number of grid points: ",i4)')
;
;  Index array:
;
kc=intarr(nc,3)
for k=0,nc-1 do begin
  dst=sqrt((xc-xc(k))^2+(yc-yc(k))^2)
  kk=where(dst lt 1.1)
  kk=kk(where(kk ne k,num))
  kc(k,0:num-1)=kk+1
endfor
;
;  Connect left and right sides (periodic boundary conditions):
;
kleft =2*indgen(jmax)+2
kright=kleft+(2*imax-1)*(2*jmax+1)
kc(kleft -1,2)=kright
kc(kright-1,2)=kleft

nt=0L  &  time=0.0
ns=0
disp_part
;
end

;-----------------------------------------------------------------
;  Read setup file:
;-----------------------------------------------------------------
pro subr_setup,error
@hexa.blk
common ran,seed
forward_function get_param
;
error=0
;
error=get_param(1,'vsetup',vsetup,/integer)
if error then return
print,string(vsetup,   $
      format='("Version of setup file:      vsetup=",i2)')


if vsetup eq 2 then begin
      ;
      error=get_param(1,'nrelax',nrelax,/integer)
      if error then return
      print,string(nrelax,   $
	    format='("Initial relaxation steps:   nrelax=",i5)')
      ;
      error=get_param(1,'nmajor',nmajor,/integer)
      if error then return
      print,string(nmajor,   $
	    format='("Major time steps:           nmajor=",i5)')
      ;
      error=get_param(1,'nminor',nminor,/integer)
      if error then return
      print,string(nminor,   $
	    format='("Minor time steps:           nminor=",i5)')
      ;
      error=get_param(1,'nrep',nrep,/integer)
      if error then return
      print,string(nrep,   $
	    format='("Repeat time steps:          nrep  =",i5)')
      ;
      error=get_param(1,'nsplit',nsplit,/float)
      if error then return
      print,string(nsplit,   $
	    format='("Splittings per major step:  nsplit=",f5.2)')
      ;
      error=get_param(1,'nbpole',nbpole,/float)
      if error then return
      print,string(nbpole,   $
	    format='("Bipoles per major step:     nbpole=",f5.2)')
      ;
      error=get_param(1,'orient',orient,/integer)
      if error then return
      print,string(orient,   $
	    format='("Bipole orientation:         orient=",i5)')
      ;
      error=get_param(1,'imode',imode,/integer)
      if error then return
      print,string(imode,   $
	    format='("Insertion mode:             imode =",i5)')
      ;
      error=get_param(1,'nsinit',nsinit,/integer)
      if error then return
      print,string(nsinit,   $
	    format='("Initial number of sources:  nsinit=",i5)')
      ;
      error=get_param(1,'seed',seed,/long)
      if error then return
      print,string(seed,   $
	    format='("Seed for random generator:  seed  =",i7)')
      ;
      error=get_param(1,'rad',rad,/float)
      if error then return
      print,string(rad,   $
	    format='("Source radius (1/e-width):  rad   =",f5.2)')
      ;
      error=get_param(1,'nfil',nfil,/integer)
      if error then return
      print,string(nfil,   $
	    format='("Insert flux rope (1=yes):   nfil  =",i2)')


endif else begin
;print,'SETUP 3'
nrelax=0
error=get_param(1,'nmajor',nmajor,/integer)
      if error then return
      print,string(nmajor,   $
	    format='("Major time steps:           nmajor=",i5)')


      nminor=1

      error=get_param(1,'nstrt',nstrt,/integer)
      if error then return
      print,string(nstrt,   $
	    format='("Start frame:           nstrt=",i5)')

      nrep=500

	    error=get_param(1,'nend',nend,/integer)
      if error then return
      print,string(nend,   $
	    format='("Ending frame:           nend=",i5)')


      nfil=0


      openr,10,dir+'param1'
      readf,10,dum
      readf,10,dum,dum,dum
      readf,10,dum
      readf,10,time_s
      close,10


      dtime=60./time_s

endelse

if nfil gt 0 then begin
;
  error=get_param(1,'ncvt',ncvt,/integer)
  if error then return
  print,string(ncvt,   $
        format='("Cavity height:           ncvt  =",i2)')
;
  error=get_param(1,'ncvw',ncvw,/integer)
  if error then return
  print,string(ncvw,   $
        format='("Cavity width:            ncvw  =",i2)')
;
  error=get_param(1,'nins',nins,/integer)
  if error then return
  print,string(nins,   $
        format='("Insertion height:        nins  =",i2)')
;
  error=get_param(1,'nwdt',nwdt,/integer)
  if error then return
  print,string(nwdt,   $
        format='("Insertion width:         nwdt  =",i2)')
;
  error=get_param(1,'faxi',faxi,/float)
  if error then return
  print,string(faxi,   $
        format='("Axial flux:              faxi  =",f7.3)')
;
  error=get_param(1,'fpol',fpol,/float)
  if error then return
  print,string(fpol,   $
        format='("Poloidal flux:           fpol  =",f7.3)')
;
endif
;
if vsetup eq 2 then begin
  error=get_param(1,'bxcon',bxcon,/float)
  if error then return
  print,string(bxcon,   $
        format='("Constant Bx:                bxcon =",f5.1," G")')
  ;
dtminor=1.0/float(nminor)    ; time step for saving (AX,AY)-values
dtime=dtminor/nrep           ; time step for magneto-friction code
endif else bxcon=0.0

end

;------------------------------------------------------------------
;  Get parameter value:
;------------------------------------------------------------------
function get_param,unit,char,var,float=float,integer=integer,long=long
;
text=''
readf,unit,text
if strlen(text) eq 0 then begin
  print,'Error in get_param: Empty string'
  return,1
endif
n=strpos(text,'=')
if n eq -1 then begin
  print,'Error in get_param: Missing equal sign'
  return,1
endif
name=strcompress(gettok(text,'='),/remove_all)
if strlen(name) eq 0 then begin
  print,'Error in get_param: No name before equal sign'
  return,1
endif
if name ne char then begin
  print,'Error in get_param: Incorrect parameter name'
  print,'Expected ',char,', found ',name
  return,1
endif
if strcompress(text,/remove_all) eq '' then begin
  print,'Error in get_param: Missing parameter value'
  return,1
endif
var=0.0
if keyword_set(float) then var=0.0
if keyword_set(integer) then var=0
if keyword_set(long) then var=0L
reads,text,var
return,0
end

;-----------------------------------------------------------------
;  Display particles on grid:
;-----------------------------------------------------------------
pro disp_part,nwin=nwin,nogrid=nogrid
@hexa.blk
;
title=string(time,ns,format='("TIME=",f5.1,", NS=",i4)')
;
if n_elements(nwin) eq 0 then nwin=0
device,window_state=flag
if flag(nwin) then wset,nwin else window,nwin,xs=500,ys=500
;
;  Display grid:
;
plot,[0.,xmax],[0.,ymax],/nodata,title=title,  $
     xrange=[0.,xmax],xstyle=1,xtitle='X',    $
     yrange=[0.,ymax],ystyle=1,ytitle='Y'
;
if not keyword_set(nogrid) then begin
  for kc1=1,nc-1 do begin
    k=kc1-1
    for i=0,2 do begin
      if kc(k,i) gt kc1 then begin
        xc2=xc(kc(k,i)-1)
        yc2=yc(kc(k,i)-1)
        dst=sqrt((xc2-xc(k))^2+(yc2-yc(k))^2)
        if dst lt 1.1 then oplot,[xc(k),xc2],[yc(k),yc2]
      endif
    endfor
  endfor
endif
;
;  Display particle positions:
;
if ns gt 0 then begin
  ind=where(fs gt 0.)
  oplot,xs(ind),ys(ind),psym=6,thick=3,color=255
  ind=where(fs le 0.)
  oplot,xs(ind),ys(ind),psym=4,thick=3,color=255*256L
endif
;
wait,0.1
;
end

;-----------------------------------------------------------------
;  Display flux density:
;-----------------------------------------------------------------
pro disp_flux,nwin=nwin
@hexa.blk
;
if n_elements(nwin) eq 0 then nwin=2
device,window_state=flag
if flag(nwin) then wset,nwin else window,nwin,xs=500,ys=500
;
nbin=10
ydel=ymax/nbin
yy=ydel*(findgen(nbin)+0.5)
fpos=fltarr(nbin)
fneg=fltarr(nbin)
for n=0,nbin-1 do begin
  y1=n*ydel
  y2=y1+ydel
  ind=where(ys ge y1 and ys lt y2 and fs gt 0,count)
  if count gt 0 then fpos(n)=total(fs(ind))
  ind=where(ys ge y1 and ys lt y2 and fs lt 0,count)
  if count gt 0 then fneg(n)=total(fs(ind))
endfor
;
;fmax=max([fpos,abs(fneg)])>150.
fmax=300.
plot,yy,fpos+fneg,title=title,  $
     xrange=[0.,ymax],xstyle=1,xtitle='Y',  $
     yrange=[-fmax,fmax],ytitle='FLUX'
;
oplot,yy,fpos,psym=10,color=255
oplot,yy,fneg,psym=10,color=255*256L
end

;-----------------------------------------------------------------
;  Compute magnetic field and vector potential at z=0:
;-----------------------------------------------------------------
pro magn,nwin=nwin
@hexa.blk
@hexa.col
;
if n_elements(nwin) eq 0 then nwin=4
device,window_state=flag
if flag(nwin) then wset,nwin else window,nwin,xs=500,ys=500
if ns eq 0 then begin
  print,'MAGN: no magnetic elements (NS=0)'
  erase
  return
endif
;
print,string(time,format='("TIME=",f5.1,", magnetic field")')
;
bz0=fltarr(nx,ny)
for n=0,ns-1 do begin
  x1=xs(n)-3.*rad  &  i1=fix(x1/delx-0.5)  &  x1=delx*(i1+0.5)
  x2=xs(n)+3.*rad  &  i2=fix(x2/delx-0.5)+1
  mx=i2-i1+1
  y1=ys(n)-3.*rad  &  j1= fix(y1/dely-0.5)>0  &  y1=dely*(j1+0.5)
  y2=ys(n)+3.*rad  &  j2=(fix(y2/dely-0.5)+1)<(ny-1)
  my=j2-j1+1
  dum=exp(-((x1+delx*findgen(mx)-xs(n))/rad)^2)  $
     #exp(-((y1+dely*findgen(my)-ys(n))/rad)^2)
  const=fs(n)/(delx*dely*total(dum))
  for m=0,mx-1 do begin
    i=i1+m
    if i lt  0 then i=i+nx
    if i ge nx then i=i-nx
    bz0(i,j1:j2)=bz0(i,j1:j2)+const*reform(dum(m,*))
  endfor
endfor
ny2=2*ny
bzd(0:nx-1, 0:ny -1)=bz0
bzd(0:nx-1,ny:ny2-1)=reverse(bz0,2)
h=float(fft(fft(bzd,-1)/qq,1))
aax0(0:nx-1,1:ny-1)= (h(0:nx-1,1:ny-1)-h(0:nx-1,0:ny-2))/dely
aay0(1:nx-1,0:ny-1)=-(h(1:nx-1,0:ny-1)-h(0:nx-2,0:ny-1))/delx
aay0(     0,0:ny-1)=-(h(     0,0:ny-1)-h(  nx-1,0:ny-1))/delx
aay0(    nx,0:ny-1)=aay0(0,0:ny-1)
aay0=aay0+(fltarr(nx+1)+1.0)#daay0     ; add contribution of the flux rope
if bxcon ne 0. then aay0=aay0+bxcon*delz*float(nz)
;
;  Write (AAX0,AAY0) on "evolution" file:
;
writeu,1,aax0
writeu,1,aay0
;
;  Display magnetogram:
;
x=delx*(findgen(nx)+0.5)
y=dely*(findgen(ny)+0.5)
bz=(aay0(1:nx,*)-aay0(0:nx-1,*))/delx-(aax0(*,1:ny)-aax0(*,0:ny-1))/dely
bmax=max(abs(bz))
title=string(time,ns,bmax,format='("TIME=",f5.1,", NS=",i4,", BMAX=",f7.2)')
bmax_plt=bmax   ; or fixed value
;
if bmax_plt ge 0.005 then begin
  nlev=5
  blev=bmax_plt*(findgen(nlev)+0.5)/nlev
  blev=[-reverse(blev),blev]
  color=[lonarr(nlev)+green,lonarr(nlev)+red]
  contour,bz,x,y,title=title,levels=blev,c_color=color,      $
          xrange=[0.,xmax],xstyle=1,xtitle='X',  $
          yrange=[0.,ymax],ystyle=1,ytitle='Y'
endif else begin
  plot,[0.,xmax],[0.,ymax],/nodata,title=title,  $
          xrange=[0.,xmax],xstyle=1,xtitle='X',  $
          yrange=[0.,ymax],ystyle=1,ytitle='Y'
endelse
;
;  Diffuse contribution from the flux rope (DAAY0) for next minor time step:
;
eta=0.05*delx^2
for n=1,nrep do begin
  diva0(1:ny-1)=(daay0(1:ny-1)-daay0(0:ny-2))/dely
  daay0=daay0+eta*(diva0(1:ny)-diva0(0:ny-1))/dely
endfor
end

;------------------------------------------------------------------
;  Compute magnetic field averaged at cell corners:
;------------------------------------------------------------------
pro prog_field
@hexa.blk
;
;  This routine calculates the magnetic field B, current density C
;  (C=curl B), and magneto-frictional velocity V (= C x B / B**2) from
;  the vector potential A.
;
;  Initialize temporary arrays for magnetic fields at cell faces
;  (BBX,BBY,BBZ), currents at cell ribs (CCX,CCY,CCZ), and values
;  at corners:
;
; bbx=fltarr(nx+1,ny+2,nz+2)
; bby=fltarr(nx+2,ny+1,nz+2)
; bbz=fltarr(nx+2,ny+2,nz+1)
;ccx=fltarr(nx+2,ny+1,nz+1)
;ccy=fltarr(nx+1,ny+2,nz+1)
;ccz=fltarr(nx+1,ny+1,nz+2)
;
;  dimension(nx+1,ny+1,nz+1) :: bx,by,bz,bb,bbm,cx,cy,cz,ch,vx,vy,vz
;
;  Compute BBX, BBY, BBZ on cell faces:
;
; bbz(1:nx,1:ny,0:nz)=(aay(1:nx,0:ny-1,0:nz)-aay(0:nx-1,0:ny-1,0:nz))/delx  $
;                    -(aax(0:nx-1,1:ny,0:nz)-aax(0:nx-1,0:ny-1,0:nz))/dely
; bbz(   0,1:ny,0:nz)=bbz(nx,1:ny,0:nz)
; bbz(nx+1,1:ny,0:nz)=bbz( 1,1:ny,0:nz)
; bbz(0:nx+1,   0,0:nz)=bbz(0:nx+1, 1,0:nz)
; bbz(0:nx+1,ny+1,0:nz)=bbz(0:nx+1,ny,0:nz)
; ;
; bbx(0:nx,1:ny,1:nz)=(aaz(0:nx,1:ny,0:nz-1)-aaz(0:nx,0:ny-1,0:nz-1))/dely  $
;                    -(aay(0:nx,0:ny-1,1:nz)-aay(0:nx,0:ny-1,0:nz-1))/delz
; bbx(0:nx,   0,1:nz)=bbx(0:nx, 1,1:nz)
; bbx(0:nx,ny+1,1:nz)=bbx(0:nx,ny,1:nz)
; bbx(0:nx,0:ny+1,  0)=bbx(0:nx,0:ny+1,   1)  $
;          -delz/delx*(bbz(1:nx+1,0:ny+1,0)-bbz(0:nx,0:ny+1,0))
; bbx(0:nx,0:ny+1,nz+1)=bbx(0:nx,0:ny+1,nz)
; ;
; bby(1:nx,0:ny,1:nz)=(aax(0:nx-1,0:ny,1:nz)-aax(0:nx-1,0:ny,0:nz-1))/delz  $
;                    -(aaz(1:nx,0:ny,0:nz-1)-aaz(0:nx-1,0:ny,0:nz-1))/delx
; bby(   0,0:ny,1:nz)=bby(nx,0:ny,1:nz)
; bby(nx+1,0:ny,1:nz)=bby(   1,0:ny,1:nz)
; bby(0:nx+1,0:ny,   0)=bby(0:nx+1,0:ny,   1)  $
;          -delz/dely*(bbz(0:nx+1,1:ny+1,0)-bbz(0:nx+1,0:ny,0))
; bby(0:nx+1,0:ny,nz+1)=bby(0:nx+1,0:ny,nz)

;
;  Compute CCX, CCY, CCZ at edges:
;
ccz=(bby(1:nx+1,0:ny  ,0:nz+1)-bby(0:nx,0:ny,0:nz+1))/delx  $
   -(bbx(0:nx  ,1:ny+1,0:nz+1)-bbx(0:nx,0:ny,0:nz+1))/dely
ccx=(bbz(0:nx+1,1:ny+1,0:nz  )-bbz(0:nx+1,0:ny,0:nz))/dely  $
   -(bby(0:nx+1,0:ny  ,1:nz+1)-bby(0:nx+1,0:ny,0:nz))/delz
ccy=(bbx(0:nx  ,0:ny+1,1:nz+1)-bbx(0:nx,0:ny+1,0:nz))/delz  $
    -(bbz(1:nx+1,0:ny+1,0:nz  )-bbz(0:nx,0:ny+1,0:nz))/delx
;
;  Magnetic field at cell corners:
;
bx1=0.25*(bbx(0:nx,0:ny,0:nz  )+bbx(0:nx,1:ny+1,0:nz  )  $
         +bbx(0:nx,0:ny,1:nz+1)+bbx(0:nx,1:ny+1,1:nz+1))
by1=0.25*(bby(0:nx,0:ny,0:nz  )+bby(1:nx+1,0:ny,0:nz  )  $
         +bby(0:nx,0:ny,1:nz+1)+bby(1:nx+1,0:ny,1:nz+1))
bz1=0.25*(bbz(0:nx,0:ny  ,0:nz)+bbz(1:nx+1,0:ny  ,0:nz)  $
         +bbz(0:nx,1:ny+1,0:nz)+bbz(1:nx+1,1:ny+1,0:nz))
bb1=bx1*bx1+by1*by1+bz1*bz1

;
;  Field strength with minimum value for magneto-friction:
;
for k=0,nz do begin
  bb1max=1.e-4*max(bb1(*,*,k))
  bb1(*,*,k)=bb1(*,*,k)>bb1max
endfor
;
;  Current density at cell corners:
;
cx1=0.5*(ccx(0:nx,0:ny,0:nz)+ccx(1:nx+1,0:ny,0:nz))
cy1=0.5*(ccy(0:nx,0:ny,0:nz)+ccy(0:nx,1:ny+1,0:nz))
cz1=0.5*(ccz(0:nx,0:ny,0:nz)+ccz(0:nx,0:ny,1:nz+1))
alpha1=(bx1*cx1+by1*cy1+bz1*cz1)/bb1
;
;  Lorentz force and magneto-frictional velocity:
;
ffx1=cy1*bz1-cz1*by1
ffy1=cz1*bx1-cx1*bz1
ffz1=cx1*by1-cy1*bx1
frc=0.1*delx^2/dtime

vx1=frc*ffx1/bb1
vy1=frc*ffy1/bb1
vz1=frc*ffz1/bb1

;
;  Magnitude of gradient(alpha):
;
ccx(1:nx,*,*)=(alpha1(1:nx,*,*)-alpha1(0:nx-1,*,*))/delx
ccx(   0,*,*)=ccx(nx,*,*)
ccx(nx+1,*,*)=ccx( 1,*,*)
ccy(*,1:ny,*)=(alpha1(*,1:ny,*)-alpha1(*,0:ny-1,*))/dely
ccy(*,   0,*)=-ccy(*, 1,*)
ccy(*,ny+1,*)=-ccy(*,ny,*)
ccz(*,*,1:nz)=(alpha1(*,*,1:nz)-alpha1(*,*,0:nz-1))/delz
ccz(*,*,   0)=-ccz(*,*, 1)
ccz(*,*,nz+1)=-ccz(*,*,nz)
grad1=sqrt((0.5*(ccx(1:nx+1,*,*)+ccx(0:nx,*,*)))^2  $
          +(0.5*(ccy(*,1:ny+1,*)+ccy(*,0:ny,*)))^2  $
          +(0.5*(ccz(*,*,1:nz+1)+ccz(*,*,0:nz)))^2)

;
;  Hyperdiffusion flux and its divergence:
;
eta4=0.01
ccx(1:nx,*,*)=0.5*eta4*(bb1(1:nx,*,*)+bb1(0:nx-1,*,*))  $
                 *(alpha1(1:nx,*,*)-alpha1(0:nx-1,*,*))/delx
ccx(   0,*,*)=ccx(nx,*,*)
ccx(nx+1,*,*)=ccx( 1,*,*)
ccy(*,1:ny,*)=0.5*eta4*(bb1(*,1:ny,*)+bb1(*,0:ny-1,*))  $
                 *(alpha1(*,1:ny,*)-alpha1(*,0:ny-1,*))/dely
ccy(*,   0,*)=-ccy(*, 1,*)
ccy(*,ny+1,*)=-ccy(*,ny,*)
ccz(*,*,1:nz)=0.5*eta4*(bb1(*,*,1:nz)+bb1(*,*,0:nz-1))  $
                 *(alpha1(*,*,1:nz)-alpha1(*,*,0:nz-1))/delz
ccz(*,*,   0)=-ccz(*,*, 1)
ccz(*,*,nz+1)=-ccz(*,*,nz)
div4=(ccx(1:nx+1,*,*)-ccx(0:nx,*,*))/delx  $
    +(ccy(*,1:ny+1,*)-ccy(*,0:ny,*))/dely  $
    +(ccz(*,*,1:nz+1)-ccz(*,*,0:nz))/delz
;
end

;------------------------------------------------------------------
;  Save vector potential on file:
;------------------------------------------------------------------
pro prog_save_model
@hexa.blk
;
;  Write file:
;
;  On SunSparc, write file in IEEE format (big endian, no extension).
;  On Intel PC (Mac), write file in native format (little endian,
;  extension "p" in file name):
;
fullname=string(filename,nt,format='(a,"_",i5.5)')
if !version.arch eq 'i386' or !version.arch eq 'x86' then begin
  extension='p'
  fullname=fullname+extension
endif
file=dir+fullname
print,'Save '+file
openw,1,file,/f77_unform
opt=1L
print,string(opt,format='("File type: opt=",i2)')
writeu,1,opt
writeu,1,aax
writeu,1,aay
writeu,1,aaz
close,1
print,'-----'
end

;------------------------------------------------------------------
;  Restore vector potential from file:
;------------------------------------------------------------------
pro prog_restore_cur
@hexa.blk
;
if filename eq '' then begin
  print,'No file selected'
  prog_restore
endif else begin
  prog_restore,string(filename,nt,format='(a,"_",i5.5)')+extension
endelse
end

;------------------------------------------------------------------
;  Restore vector potential from file:
;------------------------------------------------------------------
pro prog_restore,fullname
@hexa.blk
;
;  File names: root_nnnnn or root_nnnnnp.
;
;  Byte swapping:
;  1) Files written by Fortran on a PC (=little endian) should have
;     an extension 'p'
;  2) Files written by Fortran on a Sun (=big endian) should have
;     no extension.
;  3) IDL-written files are always in IEEE format (big endian, no extension).
;  4) other workstations (SGI, DEC) not supported.
;
;  Read file:
;
if n_elements(fullname) eq 0 then begin
  fullname=''
  if n_elements(filename) gt 0 then begin
    read,prompt='Give full name (root_##### or root_#####p) or number only:  ',$
         fullname
  endif else begin
    read,prompt='Give full name (root_##### or root_#####p):  ',fullname
  endelse
  if fullname eq '' then begin
    print,'-----'
    return
  endif
  fullname=strcompress(fullname,/remove_all)
  if n_elements(filename) gt 0 then begin
    b=byte(strmid(fullname,0,1))
    if b ge 48 and b le 57 then begin
      nt=long(fullname)
      fullname=string(filename,nt,format='(a,"_",i5.5)')+extension
    endif
  endif
endif
k=strlen(fullname)-1
if strmid(fullname,k,1) eq 'p' then begin
  extension='p'
  print,'File written by PC Fortran (extension=p)'
endif else begin
  extension=''
  print,'File in IEEE format (no extension)'
endelse
file=fullname
filename=gettok(file,'_')
nt=long(gettok(file,'p'))
print,'filename=',filename
print,'nt=',nt
file=string(filename,nt,format='(a,"_",i5.5)')+extension
if file ne fullname then begin
  print,'Error: syntax error in full name'
  print,'       file=',file,'  fullname=',fullname
  print,'-----'
  return
endif
close,1            ; make sure unit 1 is closed
;
;  Read model parameters:
;
file=dir+filename+'_setup'
openr,1,file,error=err
if err ne 0 then begin
  print,'Error: Setup file not found'
  print,'-----'
  close,1
  return
endif
print,'Read setup file: ',file
subr_setup
close,1
;
;  Read 3D data file.
;
;  Swap bytes if file was written by Fortran on a PC (extenstion='p')
;  and current program is running on a Sun workstation (=big endian).
;  Also swap bytes if file was not written by Fortran on a PC
;  (extension='') and this program is running on a PC (=little endian).
;
file=dir+fullname
if extension eq 'p' then begin
  openr,1,file,error=err,/f77_unform,/swap_if_big_endian
endif else begin
  openr,1,file,error=err,/f77_unform,/swap_if_little_endian
endelse
if err ne 0 then begin
  close,1
  print,'File '+file+' not found'
  print,'-----'
  filename=''
  nt=0L
  return
endif
opt=0L       ; read OPT parameter from file (new standard)
readu,1,opt
if opt ne 1 then begin
  close,1
  print,'Unknown format, opt=',opt
  print,'-----'
  return
endif
print,string(opt,format='("File type: opt=",i2)')
print,'Restore '+file
aax=fltarr(nx,ny+1,nz+1)
aay=fltarr(nx+1,ny,nz+1)
aaz=fltarr(nx+1,ny+1,nz)
; apkp - s
bbx = DBLARR(nx+1, ny+2, nz+2)
bby = DBLARR(nx+1, ny+2, nz+2)
bbz = DBLARR(nx+1, ny+2, nz+2)
READU,1,bbx
READU,1,bby
READU,1,bbz
; apkp - e
close,1
prog_field
end

;-------------------------------------------------------------------------
;  Create and Realize 1D Display Window:
;-------------------------------------------------------------------------
pro PROG_PLOT1D
@hexa.blk
@hexa.col
common wid,wid_main
;
;  Setup W1D structure:
;
w1d={w1d, draw: 0L, ptyp: 0L,   ds: 0L,   zh: 0L,  x1h: 0L,  y1h: 0L, $
           x2h: 0L,  y2h: 0L,   xv: 0L,   yv: 0L,  z1v: 0L,  z2v: 0L, $
          var1: 0L, cmp1: 0L, var2: 0L, cmp2: 0L, var3: 0L, cmp3: 0L, $
          ylab: 0L, nwin: 0L}
;
;  Create plot1d structure:
;
if n_elements(plot1d) eq 0 then begin
  xpos=0.5*(xmin+xmax)
  ypos=0.5*(ymin+ymax)
  plot1d={PLOT1D, upd: 0, ptyp: 0, zh: 0., x1h: xmin,  $
                  y1h: ypos, x2h: xmax, y2h: ypos, $
                  xv: xpos, yv: ypos, z1v: 0., z2v: 10.}
endif
;
;  Create base widget:
;
base=widget_base(title='1D Display',group_leader=wid_main(0),  $
                 xoffset=300,yoffset=0,/column,uvalue=w1d)
;
;  First row contains 2D image display:
;
row1=widget_base(base,/row)
w1d.draw=widget_draw(row1,/button_events,uvalue='DRAW',   $
                           x_scroll_size=800,xsize=820,   $
                           y_scroll_size=700,ysize=720)
;
;  Second row contains image control widgets:
;
row2=widget_base(base,/row)
;
;  1) First group contains [HOR,VER] droplist, VIEW and QUIT buttons:
;
acol=widget_base(row2,/column)
w1d.ptyp=cw_bgroup(acol,['Hor','Vert'],column=2,/exclusive,/no_release,  $
                   set_value=ptyp,uvalue='PTYP')
w1d.ds=cw_field(acol,title='DS=',value=0.05,/floating,xsize=6)
item=widget_button(acol,value='VIEW',uvalue='VIEW')
item=widget_button(acol,value='View in Window 0',uvalue='VIEW0')
item=widget_button(acol,value='EPS File',uvalue='EPSF')
item=widget_button(acol,value='QUIT',uvalue='QUIT')
;
;  2) Second group contains position info:
;
acol=widget_base(row2,/column)
item=widget_label(acol,value='Horizontal:')
w1d.zh=cw_field(acol,title='Z=' ,value=plot1d.zh,/float,xsize=5)
grp=widget_base(acol,/row)
w1d.x1h=cw_field(grp,title='X1=',value=plot1d.x1h,/float,xsize=5)
w1d.y1h=cw_field(grp,title='Y1=',value=plot1d.y1h,/float,xsize=5)
grp=widget_base(acol,/row)
w1d.x2h=cw_field(grp,title='X2=',value=plot1d.x2h,/float,xsize=5)
w1d.y2h=cw_field(grp,title='Y2=',value=plot1d.y2h,/float,xsize=5)
item=widget_label(acol,value='Vertical:')
grp=widget_base(acol,/row)
w1d.xv =cw_field(grp,title=' X=',value=plot1d.xv ,/float,xsize=5)
w1d.yv =cw_field(grp,title=' Y=',value=plot1d.yv ,/float,xsize=5)
grp=widget_base(acol,/row)
w1d.z1v=cw_field(grp,title='Z1=',value=plot1d.z1v,/float,xsize=5)
w1d.z2v=cw_field(grp,title='Z2=',value=plot1d.z2v,/float,xsize=5)
;
;  3) Third group contains info about variables to be plotted:
;
acol=widget_base(row2,/column)
grp=widget_base(acol,/row)
w1d.var1=cw_bgroup(grp,['B','C','V','F','Pr','Tn','Alp','No'],  $
                  column=8,/exclusive,/no_release,label_left='Solid: ', $
                  set_value=3,uvalue='VAR1')
w1d.cmp1=cw_bgroup(grp,['X','Y','Z','Abs'],column=4,  $
                  /exclusive,/no_release,label_left='--', $
                  set_value=0,uvalue='CMP1')
grp=widget_base(acol,/row)
w1d.var2=cw_bgroup(grp,['B','C','V','F','Pr','Tn','Alp','No'],  $
                  column=8,/exclusive,/no_release,label_left='Dotted:', $
                  set_value=4,uvalue='VAR2')
w1d.cmp2=cw_bgroup(grp,['X','Y','Z','Abs'],column=4,  $
                  /exclusive,/no_release,label_left='--', $
                  set_value=0,uvalue='CMP2')
grp=widget_base(acol,/row)
w1d.var3=cw_bgroup(grp,['B','C','V','F','Pr','Tn','Alp','No'],  $
                  column=8,/exclusive,/no_release,label_left='Dashed:', $
                  set_value=5,uvalue='VAR3')
w1d.cmp3=cw_bgroup(grp,['X','Y','Z','Abs'],column=4,  $
                  /exclusive,/no_release,label_left='--', $
                  set_value=0,uvalue='CMP3')
w1d.ylab=cw_field(acol,title='Label:',value='',/string,xsize=30)
;
;  Realize widgets and register with XMANAGER:
;
widget_control,/realize,base
widget_control,w1d.draw,get_value=nwin    ; get IDL window index
w1d.nwin=nwin
widget_control,base,set_uvalue=w1d    ; store W1D in UVALUE of base widget
xmanager,'PLOT1D',base,/no_block      ; /just_reg
end

;-------------------------------------------------------------------------
;  Process Events from 1D Display Window:
;-------------------------------------------------------------------------
pro PLOT1D_EVENT,ev
@hexa.blk
@hexa.trc
@hexa.col
;
widget_control,ev.id,get_uvalue=uvalue
widget_control,ev.top,get_uvalue=w1d
;
if uvalue eq 'QUIT' then begin
  widget_control,ev.top,/destroy
  return
endif
;
;  Show 1D Display:
;
if uvalue eq 'VIEW' or uvalue eq 'VIEW0' or uvalue eq 'EPSF' then begin
  if uvalue eq 'VIEW0' then begin
    device,window_state=flag
    if flag(0) then wset,nwin else window,0,xs=950,ys=720
  endif else wset,w1d.nwin
  if plot1d.upd then begin    ; update values in PLOT1D widget
    widget_control,w1d.ptyp,set_value=plot1d.ptyp
    widget_control,w1d.zh  ,set_value=plot1d.zh
    widget_control,w1d.x1h ,set_value=plot1d.x1h
    widget_control,w1d.y1h ,set_value=plot1d.y1h
    widget_control,w1d.x2h ,set_value=plot1d.x2h
    widget_control,w1d.y2h ,set_value=plot1d.y2h
    widget_control,w1d.xv  ,set_value=plot1d.xv
    widget_control,w1d.yv  ,set_value=plot1d.yv
    widget_control,w1d.z1v ,set_value=plot1d.z1v
    widget_control,w1d.z2v ,set_value=plot1d.z2v
    plot1d.upd=0
  endif
  widget_control,w1d.ptyp,get_value=ptyp &  plot1d.ptyp=ptyp
  widget_control,w1d.zh  ,get_value=zh   &  plot1d.zh  =zh
  widget_control,w1d.x1h ,get_value=x1h  &  plot1d.x1h =x1h
  widget_control,w1d.y1h ,get_value=y1h  &  plot1d.y1h =y1h
  widget_control,w1d.x2h ,get_value=x2h  &  plot1d.x2h =x2h
  widget_control,w1d.y2h ,get_value=y2h  &  plot1d.y2h =y2h
  widget_control,w1d.xv  ,get_value=xv   &  plot1d.xv  =xv
  widget_control,w1d.yv  ,get_value=yv   &  plot1d.yv  =yv
  widget_control,w1d.z1v ,get_value=z1v  &  plot1d.z1v =z1v
  widget_control,w1d.z2v ,get_value=z2v  &  plot1d.z2v =z2v
  widget_control,w1d.ds  ,get_value=ds
  case ptyp of
  0: begin    ; horizontal
     x1=float(x1h) &  y1=float(y1h)  &  z1=float(zh)
     x2=float(x2h) &  y2=float(y2h)  &  z2=float(zh)
     smax=sqrt((x2-x1)^2+(y2-y1)^2)
     if smax le 2.*ds then begin
         POP_UP,'Error','Insufficient distance between ends'
       return
     endif
     cth=(x2-x1)/smax
     sth=(y2-y1)/smax
     end
  1: begin    ; vertical
     x1=float(xv) &  y1=float(yv)  &  z1=float(z1v)
     x2=float(xv) &  y2=float(yv)  &  z2=float(z2v)
     smax=abs(z2-z1)
     if smax le 2.*ds then begin
         POP_UP,'Error','Insufficient distance between ends'
       return
     endif
     cth=1.0
     sth=0.0
     end
  endcase
  num=round(smax/ds)+1      ; number of points along extraction line
  ds=smax/(num-1)           ; adjust ds
;  print,ds
  frac=findgen(num)/(num-1)
  s=smax*frac
  xx=x1+(x2-x1)*frac
  yy=y1+(y2-y1)*frac
  zz=z1+(z2-z1)*frac
  val=fltarr(num,3)
  px=fltarr(num)
  py=fltarr(num)
  pz=fltarr(num)
  widget_control,w1d.var1,get_value=var1
  widget_control,w1d.cmp1,get_value=cmp1
  widget_control,w1d.var2,get_value=var2
  widget_control,w1d.cmp2,get_value=cmp2
  widget_control,w1d.var3,get_value=var3
  widget_control,w1d.cmp3,get_value=cmp3
  var=[var1,var2,var3]
  cmp=[cmp1,cmp2,cmp3]
  m=0
  for mm=0,2 do begin   ; loop over line styles
    if var(mm) lt 6 then begin
      par=var(mm)<3
      for n=0,num-1 do begin
        tr_magn,xx(n),yy(n),zz(n),dd,bx,by,bz,par=par
        px(n)= cth*bx+sth*by  ; for horizontal, px is parallel component
        py(n)=-sth*bx+cth*by
        pz(n)=bz
      endfor
      if var(mm) lt 4 then begin     ; MAG, CUR, VEL or FOR
        if cmp(mm) le 2 then begin
          case cmp(mm) of
          0: val(*,m)=px
          1: val(*,m)=py
          2: val(*,m)=pz
          endcase
        endif else val(*,m)=sqrt(px^2+py^2+pz^2)
      endif else begin   ; parallel pressure gradient or tension force
        pp=fltarr(num)
        hh=fltarr(num)
        for n=0,num-1 do begin
          tr_magn,xx(n),yy(n),zz(n),dd,bx,by,bz,par=0
          pp(n)=0.5*(bx^2+by^2+bz^2)
        endfor
        dp=-(pp(1:num-1)-pp(0:num-2))/ds
        val(*,m)=[dp(0),0.5*(dp(0:num-3)+dp(1:num-2)),dp(num-2)]
        if var(mm) eq 5 then begin
          if ptyp eq 0 then val(*,m)=px-val(*,m) else val(*,m)=pz-val(*,m)
        endif
      endelse
      m=m+1
    endif
    if var(mm) eq 6 then begin  ; alpha
      for n=0,num-1 do begin
        tr_magn,xx(n),yy(n),zz(n),dd,alpha,dum,dum,par=4
        val(n,m)=alpha
      endfor
      m=m+1
    endif
  endfor
  nline=m
  if nline eq 0 then begin
    print,'No lines selected'
    return
  endif
  if uvalue eq 'EPSF' then begin
    set_plot,'PS'
    epsfile=dir+'idl.eps'
    device,filename=epsfile,/encapsulated,xsize=18.,ysize=16.
  endif
  widget_control,w1d.ylab,get_value=ytitle
  if ytitle eq '' then ytitle='VALUE'
  vmin=min(val(*,0:nline-1),max=vmax)
  nexp=round(alog10(max([vmax,-vmin])))
  vexp=10.^nexp
; ytitle=ytitle+string(nexp,format='(" [x 10!U",i2,"!D]")')
  ytitle=ytitle+string(nexp,format='(" [x 10^",i2,"]")')
  plot,s,val(*,0)/vexp,linestyle=0,charsize=2.0,  $
       xr=[0.,smax],xstyle=1,xtitle='POSITION',  $
       yr=[vmin,vmax]/vexp,ytitle=ytitle
  if nline gt 1 then begin
    for m=1,nline-1 do oplot,s,val(*,m)/vexp,linestyle=m
  endif
  if uvalue eq 'EPSF' then begin
    device,/close
    print,'Plot saved on '+epsfile
    set_plot,'X'
    savfile='save.dat'
    save,filename=dir+savfile,s,val
    print,'Data saved on '+savfile
  endif
endif
;
end

;-------------------------------------------------------------------------
;  Create and Realize 2D Display Window:
;-------------------------------------------------------------------------
pro PROG_PLOT2D
@hexa.blk
@hexa.col
common wid,wid_main
;
;  Setup W2D structure:
;
w2d={w2d, draw: 0L, ptyp: 0L, xpos: 0L, xmin: 0L, xmax: 0L, $
          ypos: 0L, ymin: 0L, ymax: 0L, zpos: 0L, zmin: 0L, zmax: 0L, $
          notl: 0L, grey: 0L, cmax: 0L, nlev: 0L, cont: 0L, abso: 0L, $
          vect: 0L, vmax: 0L, null: 0L, xtra: 0L, flid: 0L, flin: 0L, $
          smin: 0L, smax: 0L, nwin: 0L, thresh: 0L}
;
;  Create base widget:
;
base=widget_base(title='2D Display',group_leader=wid_main(0),  $
                 xoffset=300,yoffset=0,/column,uvalue=w2d)
;
;  First row contains 2D image display:
;
row1=widget_base(base,/row)
w2d.draw=widget_draw(row1,/button_events,uvalue='DRAW',   $
                           x_scroll_size=800,xsize=820,   $
                           y_scroll_size=700,ysize=720)
;
;  Second row contains image control widgets:
;
row2=widget_base(base,/row)
;
;  1) First group contains [HOR,VERT] droplist, VIEW and QUIT buttons:
;
acol=widget_base(row2,/column)
w2d.ptyp=cw_bgroup(acol,['XY','XZ','YZ','SZ'],column=4,/exclusive,  $
                    /no_release,set_value=0,uvalue='PTYP')
grp=widget_base(acol,/row)
item=widget_button(grp,value='Zoom',uvalue='ZOOM')
item=widget_button(grp,value='Fpol',uvalue='FPOL')
item=widget_button(acol,value='VIEW',uvalue='VIEW')
item=widget_button(acol,value='F-L Reset (ID=0)',uvalue='RESET')
item=widget_button(acol,value='View in Window 0',uvalue='VIEW0')
item=widget_button(acol,value='EPS File',uvalue='EPSF')
item=widget_button(acol,value='QUIT',uvalue='QUIT')
;
;  2) Second group contains position info:
;
acol=widget_base(row2,/column)
item=widget_label(acol,value='  POS        MIN       MAX     ')
grp=widget_base(acol,/row)
xpos=0.5*(xmin+xmax)
w2d.xpos=cw_field(grp,title='X: ',value=xpos,/float,xsize=5)
w2d.xmin=cw_field(grp,title=' '  ,value=xmin,/float,xsize=5)
w2d.xmax=cw_field(grp,title=' '  ,value=xmax,/float,xsize=5)
item=widget_button(grp,value='Rs',uvalue='XRES')
grp=widget_base(acol,/row)
ypos=0.5*(ymin+ymax)
w2d.ypos=cw_field(grp,title='Y: ',value=ypos,/float,xsize=5)
w2d.ymin=cw_field(grp,title=' '  ,value=ymin,/float,xsize=5)
w2d.ymax=cw_field(grp,title=' '  ,value=ymax,/float,xsize=5)
item=widget_button(grp,value='Rs',uvalue='YRES')
grp=widget_base(acol,/row)
w2d.zpos=cw_field(grp,title='Z: ',value= 0.,/float,xsize=5)
w2d.zmin=cw_field(grp,title=' '  ,value= 0.,/float,xsize=5)
w2d.zmax=cw_field(grp,title=' '  ,value=10.,/float,xsize=5)
grp=widget_base(acol,/row)
w2d.notl=cw_bgroup(grp,['Yes','No'],column=2,  $
                  /exclusive,/no_release,label_left='Title:', $
                  set_value=0,uvalue='NOTL')
item=widget_button(grp,value='F-L Alpha',uvalue='ALPHA')
;
;  3) Third group contains info for contours, vectors, etc.:
;
acol=widget_base(row2,/column)
grp=widget_base(acol,/row)
w2d.cont=cw_bgroup(grp,['B','C','V','F','Alp','Grd','Div'],column=7,  $
                  /exclusive,/no_release,label_left='Var:', $
                  set_value=0,uvalue='CONT')
w2d.abso=cw_bgroup(grp,['Cmp','Abs'],column=2,  $
                  /exclusive,/no_release,label_left='-', $
                  set_value=0,uvalue='ABSO')
grp=widget_base(acol,/row)
;w2d.grey=cw_bgroup(grp,['CC','BW','Grey','QSL','Con','No'],column=6,  $
w2d.grey=cw_bgroup(grp,['CC','BW','Grey','No'],column=4,  $
                  /exclusive,/no_release,label_left='',  $
                  set_value=0,uvalue='GREY')
w2d.cmax=cw_field(grp,title='Max:',value=0.,/floating,xsize=8)
w2d.nlev=cw_field(grp,title='Lev:',value=5,/integer,xsize=2)
grp=widget_base(acol,/row)
w2d.vect=cw_bgroup(grp,['MAG','CUR','VEL','FOR','None'],column=5,  $
                  /exclusive,/no_release,label_left='Vectors:',  $
                  set_value=0,uvalue='VECT')
w2d.vmax=cw_field (grp,title='Max:',value=0.,/floating,xsize=8)
grp=widget_base(acol,/row)
w2d.flid=cw_field(grp,title='F-L ID:',value=0,/integer,xsize=3)
w2d.flin=cw_field(grp,title='F-L Index:',value=-1,/integer,xsize=3)
w2d.smin=cw_field(grp,title='SMin:',value=-999.,/float,xsize=6)
w2d.smax=cw_field(grp,title='SMax:',value= 999.,/float,xsize=6)
item=widget_button(grp,value='Rs',uvalue='SRES')
grp=widget_base(acol,/row)
w2d.xtra=cw_bgroup(grp,['No','Yes'],column=2,  $
                  /exclusive,/no_release,label_left='Plot1D:', $
                  set_value=0,uvalue='XTRA')
item=widget_button(grp,value='Sel',uvalue='PSEL')

item=widget_button(grp,value='Find Flux Ropes',uvalue='frfind')
w2d.thresh=cw_field(grp,title='Max Field Strength',value=50.,/float,xsize=4)



;w2d.null=cw_bgroup(grp,['No','Yes'],column=2,  $
;                   /exclusive,/no_release,label_left='Nulls:',  $
 ;                  set_value=0,uvalue='NULL')
;item=widget_button(grp,value='Fan/Spine',uvalue='NSEL')

;  Realize widgets and register with XMANAGER:
;
widget_control,/realize,base
widget_control,w2d.draw,get_value=nwin    ; get IDL window index
w2d.nwin=nwin
widget_control,base,set_uvalue=w2d        ; store W2D in UVALUE of base widget
plot2d,nwin=w2d.nwin                      ; create default display
xmanager,'PLOT2D',base,/no_block   ; /just_reg
end

;-------------------------------------------------------------------------
;  Process Events from 2D Display Window:
;-------------------------------------------------------------------------
pro PLOT2D_EVENT,ev
@hexa.blk
@hexa.trc
@hexa.col
;
widget_control,ev.id,get_uvalue=uvalue
widget_control,ev.top,get_uvalue=w2d
;
if uvalue eq 'QUIT' then begin
  widget_control,ev.top,/destroy
  return
endif

if uvalue eq 'frfind' then begin
  widget_control,w2d.thresh,get_value=thresh
  ;print,thresh
  frfinder,threshold=thresh
  print,"To view field lines, please click 'VIEW'"
  ;uvalue='VIEW'
  return
endif


;
;  Reset XRANGE and YRANGE:
;
if uvalue eq 'XRES' then begin
  widget_control,w2d.xmin,set_value=xmin
  widget_control,w2d.xmax,set_value=xmax
  return
endif
if uvalue eq 'YRES' then begin
  widget_control,w2d.ymin,set_value=ymin
  widget_control,w2d.ymax,set_value=ymax
  return
endif
if uvalue eq 'SRES' then begin
  widget_control,w2d.smin,set_value=-999.
  widget_control,w2d.smax,set_value= 999.
  return
endif
;
;  Plot ALPHA along last-traced field line:
;
if uvalue eq 'ALPHA' then begin
  num=n_elements(xtrc)
  if num gt 0 then begin
    print,string(ztrc(kp),format='("Height of starting point: z=",f6.3)')
    z0=1.5*delz
    title=string(z0,format='("Dashed: z < ",f4.2)')
    slen=strc(num-1)-strc(0)
    smin=strc(    0)-0.05*slen
    smax=strc(num-1)+0.05*slen
    device,window_state=flag
    if flag(12) then wset,12 else window,12,xs=700,ys=600
    atrc=fltarr(num)
    for k=0,num-1 do begin
      tr_magn,xtrc(k),ytrc(k),ztrc(k),dum,cx2,cy2,cz2,par=1   ; get current
      atrc(k)=(cx2*bx2(k)+cy2*by2(k)+cz2*bz2(k))/bm2(k)^2
    endfor
    plot,strc,atrc,linestyle=2,title=title,   $
      xtitle='POSITION ALONG FIELD LINE',xstyle=1,xr=[smin,smax],  $
      ytitle='ALPHA'
    k=where(ztrc gt z0,count)
    if count ge 2 then begin
      dk=[k(1:count-1)-k(0:count-2),1L]
      i2=-1
      repeat begin
        i1=i2+1
        i2=i1
        while (dk(i2) eq 1 and i2 le (count-2)) do i2=i2+1
        kk=k(i1:i2)
;       print,'kk=',kk
        oplot,strc(kk),atrc(kk)
      endrep until (i2 eq (count-1))
    endif
;   !p.multi=0
  endif
  uvalue='VIEW'
endif
;
;  Select Zoom region:
;
if uvalue eq 'ZOOM' then begin
  if w2d.nwin ne !d.window then begin
    print,'First VIEW before selecting Zoom region'
  endif else begin
    widget_control,w2d.ptyp,get_value=ptyp
    if ptyp le 2 then text='to select Zoom region'  $
                 else text='to measure poloidal flux'
    print,'Click LEFT in PLOT2D window '+text+' (RIGHT=quit)'
    case ptyp of
    0: begin   ; XY
       cursor,x1,y1,3,/data
       if !err eq 4 then return
;       x1=round(x1)  &  y1=round(y1)
       cursor,x2,y2,3,/data
       if !err eq 4 then return
;       x2=round(x2)  &  y2=round(y2)
       x1=min([x1,x2],max=x2)
       y1=min([y1,y2],max=y2)
       widget_control,w2d.xmin,set_value=x1
       widget_control,w2d.xmax,set_value=x2
       widget_control,w2d.ymin,set_value=y1
       widget_control,w2d.ymax,set_value=y2
       end
    1: begin   ; XZ
       cursor,x1,z1,3,/data
;       x1=round(x1)  &  z1=round(z1)
       if !err eq 4 then return
       cursor,x2,z2,3,/data
       if !err eq 4 then return
;       x2=round(x2)  &  z2=round(z2)
       x1=min([x1,x2],max=x2)
       z1=min([z1,z2],max=z2)
       widget_control,w2d.xmin,set_value=x1
       widget_control,w2d.xmax,set_value=x2
       widget_control,w2d.zmin,set_value=z1
       widget_control,w2d.zmax,set_value=z2
       end
    2: begin   ; YZ
       cursor,y1,z1,3,/data
;       y1=round(y1)  &  z1=round(z1)
       if !err eq 4 then return
       cursor,y2,z2,3,/data
       if !err eq 4 then return
;       y2=round(y2)  &  z2=round(z2)
       y1=min([y1,y2],max=y2)
       z1=min([z1,z2],max=z2)
       widget_control,w2d.ymin,set_value=y1
       widget_control,w2d.ymax,set_value=y2
       widget_control,w2d.zmin,set_value=z1
       widget_control,w2d.zmax,set_value=z2
       end
    3: begin   ; SZ
       POP_UP,'Error','Cannot Zoom in SZ display'
       return
       end
    endcase
  endelse
  uvalue='VIEW'
endif
;
;  Measure poloidal flux:
;
if uvalue eq 'FPOL' then begin
  if w2d.nwin ne !d.window then begin
    POP_UP,'Error','First VIEW this window'
    return
  endif
  widget_control,w2d.ptyp,get_value=ptyp
  if ptyp le 2 then begin
    POP_UP,'Error','Poloidal Flux measurement only works in SZ display'
    return
  endif else begin
    print,'Click LEFT in PLOT2D window (RIGHT=quit)'
    widget_control,w2d.zmin,get_value=zmini
    widget_control,w2d.zmax,get_value=zmaxi
    zmini=float(zmini)+0.01
    zmaxi=float(zmaxi)-0.01
    x1h=float(plot1d.x1h)
    y1h=float(plot1d.y1h)
    x2h=float(plot1d.x2h)
    y2h=float(plot1d.y2h)
    smax=sqrt((x2h-x1h)^2+(y2h-y1h)^2)
    cth=(x2h-x1h)/smax
    sth=(y2h-y1h)/smax
    cursor,s1,z1,3,/data
    if !err eq 4 then return
    s1=s1>0.<smax
    x1=x1h+s1*cth
    y1=y1h+s1*sth
    z1=z1>zmini<zmaxi
    cursor,s2,z2,3,/data
    if !err eq 4 then return
    s2=s2>0.<smax
    x2=x1h+s2*cth
    y2=y1h+s2*sth
    z2=z2>zmini<zmaxi
    len=sqrt((s2-s1)^2+(z2-z1)^2)
    if len lt 1.0 then begin
       POP_UP,'Error','Selected interval is too short'
       return
    endif
    num=round(len)   ; number of intervals
    ds=(s2-s1)/num
    dz=(z2-z1)/num
    fpol=0.0
    for n=0,num-1 do begin
      pos=float(n)+0.5
      s=s1+pos*ds
      x=x1h+s*cth
      y=y1h+s*sth
      z=z1+pos*dz
      tr_magn,x,y,z,dum,bx2,by2,bz2,hx2,hy2,hz2
      len=sqrt((ds*hx2)^2+(dz*hz2)^2)
      fpol=fpol+(ds*hx2*bz2-dz*hz2*(cth*bx2+sth*by2))*rsun_cm
    endfor
    print,string(fpol,format='("Poloidal flux: ",e10.3," Mx/cm")')
    np0=2
    xp0=[x1,x2]
    yp0=[y1,y2]
    zp0=[z1,z2]
    widget_control,w2d.null,set_value=1   ; show as "nulls"
  end
  uvalue='VIEW'
endif
;
;  Select PLOT1D extraction line:
;
if uvalue eq 'PSEL' then begin
  if w2d.nwin ne !d.window then begin
    POP_UP,'Error','First VIEW this window before selecting PLOT1D line'
    return
  endif
  if n_elements(plot1d) eq 0 then begin
    plot1d={PLOT1D, upd: 0, ptyp: 0,  $
                    zh: 0., x1h: 0., y1h: 0., x2h: 0., y2h: 0., $
                    xv: 0., yv: 0., z1v: 0., z2v: 0.}
  endif
  widget_control,w2d.ptyp,get_value=ptyp
  if ptyp eq 0 then begin
    print,'Click LEFT in PLOT2D to select horizontal line (RIGHT=quit)'
    widget_control,w2d.zpos,get_value=z1
    cursor,x1,y1,3,/data
    if !err eq 4 then return
    x1=x1>xmin<xmax
    y1=y1>ymin<ymax
    print,string(x1,y1,z1,format='("x1=",f7.1," y1=",f7.1," z1=",f7.1)')
    cursor,x2,y2,3,/data
    if !err eq 4 then return
    x2=x2>xmin<xmax
    y2=y2>ymin<ymax
    print,string(x2,y2,z1,format='("x2=",f7.1," y2=",f7.1," z2=",f7.1)')
    plot1d.ptyp=0
    plot1d.zh =z1
    plot1d.x1h=x1
    plot1d.y1h=y1
    plot1d.x2h=x2
    plot1d.y2h=y2
  endif else begin
    print,'Click LEFT in PLOT2D to select vertical line (RIGHT=quit)'
    case ptyp of
    1: begin   ; XZ
       widget_control,w2d.ypos,get_value=yy
       y1=float(yy)
       cursor,x1,z1,3,/data
       x1=x1>xmin<xmax
       z1=z1>zmin<zmax
       if !err eq 4 then return
       print,string(x1,y1,z1,format='("x1=",f7.1," y1=",f7.1," z1=",f7.1)')
       cursor,dum,z2,3,/data
       if !err eq 4 then return
       z2=z2>zmin<zmax
       print,string(x1,y1,z2,format='("x2=",f7.1," y2=",f7.1," z2=",f7.1)')
       end
    2: begin   ; YZ
       widget_control,w2d.xpos,get_value=xx
       x1=float(xx)
       cursor,y1,z1,3,/data
       if !err eq 4 then return
       y1=y1>ymin<ymax
       z1=z1>zmin<zmax
       print,string(x1,y1,z1,format='("x1=",f7.1," y1=",f7.1," z1=",f7.1)')
       cursor,dum,z2,3,/data
       if !err eq 4 then return
       z2=z2>zmin<zmax
       print,string(x1,y1,z2,format='("x2=",f7.1," y2=",f7.1," z2=",f7.1)')
       end
    3: begin   ; SZ
       cursor,s1,z1,3,/data
       if !err eq 4 then return
       x1h=float(plot1d.x1h)
       y1h=float(plot1d.y1h)
       x2h=float(plot1d.x2h)
       y2h=float(plot1d.y2h)
       smax=sqrt((x2h-x1h)^2+(y2h-y1h)^2)
       cth=(x2h-x1h)/smax
       sth=(y2h-y1h)/smax
       x1=x1h+s1*cth
       y1=y1h+s1*sth
       x1=x1>xmin<xmax
       y1=y1>ymin<ymax
       z1=z1>zmin<zmax
       print,string(x1,y1,z1,format='("x1=",f7.1," y1=",f7.1," z1=",f7.1)')
       cursor,dum,z2,3,/data
       if !err eq 4 then return
       z2=z2>zmin<zmax
       print,string(x1,y1,z2,format='("x2=",f7.1," y2=",f7.1," z2=",f7.1)')
       end
    endcase
    z1=min([z1,z2],max=z2)
    plot1d.ptyp=1
    plot1d.xv =x1
    plot1d.yv =y1
    plot1d.z1v=z1
    plot1d.z2v=z2
  endelse
  plot1d.upd=1
  widget_control,w2d.xtra,set_value=1
  uvalue='VIEW'
endif
;
;  Select nearest null:
;
if uvalue eq 'NSEL' then begin
  if np0 eq 0 then begin
    POP_UP,'Error','No nulls selected'
    return
  endif
  if w2d.nwin ne !d.window then begin
    POP_UP,'Error','First VIEW this window before selecting null'
    return
  endif
  print,'Click LEFT in PLOT2D to select a NULL point (yellow cross)'
  widget_control,w2d.ptyp,get_value=ptyp
  case ptyp of
  0: begin
     cursor,x1,y1,3,/data
     if !err eq 4 then return
     dstmin=min(sqrt((xp0-x1)^2+(yp0-y1)^2),nn)
     end
  1: begin
     cursor,x1,z1,3,/data
     if !err eq 4 then return
     dstmin=min(sqrt((xp0-x1)^2+(zp0-z1)^2),nn)
     end
  2: begin
     cursor,y1,z1,3,/data
     if !err eq 4 then return
     dstmin=min(sqrt((yp0-y1)^2+(zp0-z1)^2),nn)
     end
  3: begin
     cursor,s1,z1,3,/data
     if !err eq 4 then return
     x1h=float(plot1d.x1h)
     y1h=float(plot1d.y1h)
     x2h=float(plot1d.x2h)
     y2h=float(plot1d.y2h)
     smax=sqrt((x2h-x1h)^2+(y2h-y1h)^2)
     cth=(x2h-x1h)/smax
     sth=(y2h-y1h)/smax
     sp0=(xp0-x1h)*cth+(yp0-y1h)*sth
     dstmin=min(sqrt((sp0-s1)^2+(zp0-z1)^2),nn)
     end
  endcase
  print,string(nn,xp0(nn),yp0(nn),zp0(nn),  $
        format='("nn=",i5,",  xp0=",f7.3,",  yp0=",f7.3,",  zp0=",f7.3)')
  prog_null_launch,nn
  uvalue='VIEW'
endif
;
;  Show 2D Display:
;
if uvalue eq 'VIEW' or uvalue eq 'VIEW0' or  $
   uvalue eq 'EPSF' or uvalue eq 'RESET' then begin
  if uvalue eq 'RESET' then subr_reset        ; reset field-line counter
  if uvalue eq 'VIEW0' then nwin=0 else nwin=w2d.nwin
  if uvalue eq 'EPSF' then epsfile=dir+'idl.eps' else epsfile=''
  widget_control,w2d.ptyp,get_value=ptyp
  widget_control,w2d.xmin,get_value=xrmin
  widget_control,w2d.xmax,get_value=xrmax
  xrange=[xrmin,xrmax]
  widget_control,w2d.ymin,get_value=yrmin
  widget_control,w2d.ymax,get_value=yrmax
  yrange=[yrmin,yrmax]
  widget_control,w2d.zmin,get_value=zrmin
  widget_control,w2d.zmax,get_value=zrmax
  zrange=[zrmin,zrmax]
  case ptyp of
  0: begin
     xz=0 & yz=0 & sz=0    ; XY
     widget_control,w2d.zpos,get_value=pos
     end
  1: begin
     xz=1 & yz=0 & sz=0    ; XZ
     widget_control,w2d.ypos,get_value=pos
     end
  2: begin
     xz=0 & yz=1 & sz=0    ; YZ
     widget_control,w2d.xpos,get_value=pos
     end
  3: begin
     xz=0 & yz=0 & sz=1    ; SZ (arbitrary vertical slice)
     if n_elements(plot1d) eq 0 then begin
       widget_control,w2d.xpos,get_value=xpos
       widget_control,w2d.ypos,get_value=ypos
       plot1d={PLOT1D, upd: 1, ptyp: 0, zh: 0.,  $
                  x1h: xrmin, y1h: ypos, x2h: xrmax, y2h: ypos, $
                  xv: xpos, yv: ypos, z1v: zrmin, z2v: zrmax}
     endif
     if plot1d.x1h eq plot1d.x2h and plot1d.y1h eq plot1d.y2h then begin
       print,'Warning: x1h and x2h changed in PLOT1D structure'
       plot1d.x1h=round((plot1d.x1h-10)>xmin)
       plot1d.x2h=round((plot1d.x1h+20)<xmax)
       plot1d.upd=1
     endif
     end
  endcase
  widget_control,w2d.notl,get_value=notitle
  widget_control,w2d.grey,get_value=greyscale
  widget_control,w2d.cmax,get_value=cm1
  if cm1 gt 0. then cm=cm1
  widget_control,w2d.nlev,get_value=nlev
  widget_control,w2d.cont,get_value=cont
  widget_control,w2d.abso,get_value=abso
  widget_control,w2d.vect,get_value=vect
  if vect eq 4 then begin & vect=0 & novect=1 & endif
  widget_control,w2d.vmax,get_value=vm1
  if vm1 gt 0. then vm=vm1
  ;widget_control,w2d.null,get_value=null
  widget_control,w2d.xtra,get_value=xtra
  widget_control,w2d.flid,get_value=flid
  widget_control,w2d.flin,get_value=flin
  widget_control,w2d.smin,get_value=smin
  widget_control,w2d.smax,get_value=smax
  plot2d,pos,xz=xz,yz=yz,sz=sz,                                 $ ; views
           xrange=xrange,yrange=yrange,zrange=zrange,           $ ; xyz range
           cont=cont,cm=cm,nlev=nlev,greyscale=greyscale,       $ ; contours
           abso=abso,vect=vect,vm=vm,novect=novect,             $ ; vectors
           null=null,xtra=xtra,flid=flid,flin=flin,smin=smin,smax=smax,   $
           nwin=nwin,notitle=notitle,epsfile=epsfile
  return
endif
;
;  Trace field lines (click in DRAW widget):
;
if uvalue eq 'DRAW' then begin
  if ev.press ne 1 then return   ; only accept Left Button events
  if w2d.nwin ne !d.window then begin
    POP_UP,'Error','First VIEW this window before tracing field lines'
    return
  endif
; print,'Event: ',ev.x,ev.y
  widget_control,w2d.flid,get_value=flid
  if flid gt 0 then begin
    POP_UP,'Error','Cannot add field lines to file, '  $
                   +string(flid,format='("FLID=",i3)')
    return
  endif
  widget_control,w2d.ptyp,get_value=ptyp
  case ptyp of
  0: begin   ; XY
       xx=float(ev.x)/!d.x_size
       yy=float(ev.y)/!d.y_size
       if xx le !x.window(0) or xx ge !x.window(1) or  $
          yy le !y.window(0) or yy ge !y.window(1) then return
       xx=(xx-!x.s(0))/!x.s(1)
       yy=(yy-!y.s(0))/!y.s(1)
       widget_control,w2d.zpos,get_value=zz
     end
  1: begin   ; XZ
       xx=float(ev.x)/!d.x_size
       zz=float(ev.y)/!d.y_size
       if xx le !x.window(0) or xx ge !x.window(1) or  $
          zz le !y.window(0) or zz ge !y.window(1) then return
       xx=(xx-!x.s(0))/!x.s(1)
       zz=(zz-!y.s(0))/!y.s(1)
       widget_control,w2d.ypos,get_value=yy
     end
  2: begin   ; YZ
       yy=float(ev.x)/!d.x_size
       zz=float(ev.y)/!d.y_size
       if yy le !x.window(0) or yy ge !x.window(1) or  $
          zz le !y.window(0) or zz ge !y.window(1) then return
       yy=(yy-!x.s(0))/!x.s(1)
       zz=(zz-!y.s(0))/!y.s(1)
       widget_control,w2d.xpos,get_value=xx
     end
  3: begin   ; SZ
       ss=float(ev.x)/!d.x_size
       zz=float(ev.y)/!d.y_size
       if ss le !x.window(0) or ss ge !x.window(1) or  $
          zz le !y.window(0) or zz ge !y.window(1) then return
       ss=(ss-!x.s(0))/!x.s(1)
       zz=(zz-!y.s(0))/!y.s(1)
       x1=float(plot1d.x1h)
       y1=float(plot1d.y1h)
       x2=float(plot1d.x2h)
       y2=float(plot1d.y2h)
       smax=sqrt((x2-x1)^2+(y2-y1)^2)
       cth=(x2-x1)/smax
       sth=(y2-y1)/smax
       xx=x1+ss*cth
       yy=y1+ss*sth
     end
  endcase
; print,string(xx,yy,zz,format='("xx=",f6.1,", yy=",f6.1,", zz=",f6.1)')
  if n_elements(np) eq 0 then subr_reset
  if np eq 0 then begin
    xp(0)=xx
    yp(0)=yy
    zp(0)=zz
  endif else begin
    xp=[xp,xx]
    yp=[yp,yy]
    zp=[zp,zz]
  endelse
  np=np+1
  widget_control,w2d.smin,get_value=smin
  widget_control,w2d.smax,get_value=smax
  trace,xx,yy,zz,smin=smin,smax=smax
  case ptyp of
  0: oplot,xtrc,ytrc,color=white,thick=2
  1: oplot,xtrc,ztrc,color=white,thick=2
  2: oplot,ytrc,ztrc,color=white,thick=2
  3: oplot,(xtrc-x1)*cth+(ytrc-y1)*sth,ztrc,color=white,thick=2
  endcase
endif
end

;------------------------------------------------------------------
;  Plot magnetic field in plane "pos":
;------------------------------------------------------------------
pro plot2d,pos,xz=xz,yz=yz,sz=sz,                            $ ; views
           xrange=xrange,yrange=yrange,zrange=zrange,        $ ; xyz range
           cont=cont,cm=cm,nlev=nlev,greyscale=greyscale, $
           abso=abso,vect=vect,vm=vm,novect=novect,               $
           flid=flid,flin=flin,smin=smin,smax=smax,nwin=nwin,null=null,     $
           xtra=xtra,notitle=notitle,epsfile=epsfile
@hexa.blk
@hexa.trc
;@hexa.qsl
@hexa.col
;
;  Setup plot device:
;
if keyword_set(epsfile) then begin
  set_plot,'PS'
  device,/encapsulated,/portrait,/color,bits=8,  $
         filename=epsfile,xsize=18.,ysize=15.8
  hexa_colour
  cc1=black
  cc2=black
  cc3=black
  cc4=red
  charsize=0.9  ; for plot title only
  th2bw=6.     ; whiteout line thickness for field lines
  th1bw=5.      ; drawing line thickness for field lines
endif else begin
  set_plot,'X'
  if n_elements(nwin) eq 0 then nwin=0
  if nwin gt 0 then wset,nwin else window,0,xsize=820,ysize=720
  hexa_colour
  cc1=white
  cc2=black
  cc3=black    ; not used
  cc4=yellow    ; magnetic nulls
  charsize=1.5  ; for plot title only
  th2bw=6.      ; whiteout line thickness for field lines
  th1bw=2.      ; drawing line thickness for field lines
endelse
ratwin=float(!d.y_size)/float(!d.x_size)
;
;  Define plot range:
;
if not keyword_set(xrange) then xrange=[xmin,xmax]
xrange(0)=xrange(0)>xmin<xmax
xrange(1)=xrange(1)>xmin<xmax
if (xrange(1)-xrange(0)) lt 2.*delx then begin
  print,'Error: invalid XRANGE'
  print,'-----'
  return
endif
if not keyword_set(yrange) then yrange=[ymin,ymax]
yrange(0)=yrange(0)>ymin<ymax
yrange(1)=yrange(1)>ymin<ymax
if (yrange(1)-yrange(0)) lt 2.*dely then begin
  print,'Error: invalid YRANGE'
  print,'-----'
  return
endif
if not keyword_set(zrange) then zrange=[zmin,zmax]
zrange(0)=zrange(0)>zmin<zmax
zrange(1)=zrange(1)>zmin<zmax
if (zrange(1)-zrange(0)) lt 2.*delz then begin
  print,'Error: invalid ZRANGE'
  print,'-----'
  return
endif
;
;  Type of vectors:
;
if n_elements(vect) eq 0 then vect=0
vect=fix(vect)>0<3
vchr=['V_MAG','V_CUR','V_VEL','V_FOR']
;
;  Type of contours:
;
if n_elements(cont) eq 0 then cont=0
cont=fix(cont)>0<6
cchr=['C_MAG','C_CUR','C_VEL','C_FOR','C_ALPHA','C_GRAD','C_DIV4']
if n_elements(abso) eq 0 then abso=0
abso=fix(abso)>0<1
if cont eq 4 then abso=0   ; display ALPHA as if a vector component
if cont eq 5 then abso=1   ; display GRAD  as if absolute value of vector
if cont eq 6 then abso=0   ; display DIV4  as if a vector component
;
;  Type of plot:
;
ptyp=0
if keyword_set(xz) then ptyp=1
if keyword_set(yz) then ptyp=2
if keyword_set(sz) and n_elements(plot1d) gt 0 then ptyp=3
;
;  Extract plot data:
;
if ptyp gt 0 then begin   ; Show vertical slices
  vchr=vchr(vect)         ; all vector  types allowed
  cchr=cchr(cont)         ; all contour types allowed
  chr='MOD'
  case ptyp of
  1: begin           ; XZ plane
     if n_elements(pos) eq 0 then pos=0.5*(ymin+ymax)
     pos=pos>ymin<ymax
     ky=round((pos-ymin)/dely)>0<ny
     pos=ymin+dely*float(ky)
     nx0=round((xrange(0)-xmin)/delx)<(nx-1)  &  xrange(0)=xmin+delx*float(nx0)
     nx1=round((xrange(1)-xmin)/delx)<nx      &  xrange(1)=xmin+delx*float(nx1)
     nz0=round((zrange(0)-zmin)/delz)<(nz-1)  &  zrange(0)=zmin+delz*float(nz0)
     nz1=round((zrange(1)-zmin)/delz)<nz      &  zrange(1)=zmin+delz*float(nz1)
     nxn=nx1-nx0+1
     nzn=nz1-nz0+1
     if cont eq 0 or vect eq 0 then begin
       bx3=reform(bx1(nx0:nx1,ky,nz0:nz1))
       by3=reform(by1(nx0:nx1,ky,nz0:nz1))
       bz3=reform(bz1(nx0:nx1,ky,nz0:nz1))
     endif
     if cont eq 1 or vect eq 1 then begin
       cx3=reform(cx1(nx0:nx1,ky,nz0:nz1))
       cy3=reform(cy1(nx0:nx1,ky,nz0:nz1))
       cz3=reform(cz1(nx0:nx1,ky,nz0:nz1))
     endif
     if cont eq 2 or vect eq 2 then begin
       vx3=reform(vx1(nx0:nx1,ky,nz0:nz1))
       vy3=reform(vy1(nx0:nx1,ky,nz0:nz1))
       vz3=reform(vz1(nx0:nx1,ky,nz0:nz1))
     endif
     if cont eq 3 or vect eq 3 then begin
       fx3=reform(ffx1(nx0:nx1,ky,nz0:nz1))
       fy3=reform(ffy1(nx0:nx1,ky,nz0:nz1))
       fz3=reform(ffz1(nx0:nx1,ky,nz0:nz1))
     endif
     if cont eq 4 then b3=reform(alpha1(nx0:nx1,ky,nz0:nz1))
     if cont eq 5 then b3=reform( grad1(nx0:nx1,ky,nz0:nz1))
     if cont eq 6 then b3=reform(  div4(nx0:nx1,ky,nz0:nz1))
     if cont le 3 then begin
       if abso then begin
         case cont of
         0: b3=sqrt(bx3^2+by3^2+bz3^2)
         1: b3=sqrt(cx3^2+cy3^2+cz3^2)
         2: b3=sqrt(vx3^2+vy3^2+vz3^2)
         3: b3=sqrt(fx3^2+fy3^2+fz3^2)
         endcase
         cchr=cchr+'_A'
       endif else begin
         case cont of
         0: b3=by3
         1: b3=cy3
         2: b3=vy3
         3: b3=fy3
         endcase
         cchr=cchr+'_Y'
       endelse
     endif
     case vect of
     0: begin & b1=bx3 & b2=bz3 & end
     1: begin & b1=cx3 & b2=cz3 & end
     2: begin & b1=vx3 & b2=vz3 & end
     3: begin & b1=fx3 & b2=fz3 & end
     endcase
     v1=xrange(0)+delx*findgen(nxn)
     v2=zrange(0)+delz*findgen(nzn)
     vchr=vchr+'_XZ'
     xtitle='X'
     ytitle='Z'
     pchr='Y='
     xrange_plt=[xrange(0)-0.5*delx,xrange(1)+0.5*delx]
     yrange_plt=[zrange(0)-0.5*delz,zrange(1)+0.5*delz]
     ratio=float(nzn)/nxn
     end
  2: begin      ; YZ plane
     if n_elements(pos) eq 0 then pos=0.5*(xmin+xmax)
     pos=pos>xmin<xmax
     kx=round((pos-xmin)/delx)>0<nx
     pos=xmin+delx*float(kx)
     ny0=round((yrange(0)-ymin)/dely)<(ny-1)  &  yrange(0)=ymin+dely*float(ny0)
     ny1=round((yrange(1)-ymin)/dely)<ny      &  yrange(1)=ymin+dely*float(ny1)
     nz0=round((zrange(0)-zmin)/delz)<(nz-1)  &  zrange(0)=zmin+delz*float(nz0)
     nz1=round((zrange(1)-zmin)/delz)<nz      &  zrange(1)=zmin+delz*float(nz1)
     nyn=ny1-ny0+1
     nzn=nz1-nz0+1
     if cont eq 0 or vect eq 0 then begin
       bx3=reform(bx1(kx,ny0:ny1,nz0:nz1))
       by3=reform(by1(kx,ny0:ny1,nz0:nz1))
       bz3=reform(bz1(kx,ny0:ny1,nz0:nz1))
     endif
     if cont eq 1 or vect eq 1 then begin
       cx3=reform(cx1(kx,ny0:ny1,nz0:nz1))
       cy3=reform(cy1(kx,ny0:ny1,nz0:nz1))
       cz3=reform(cz1(kx,ny0:ny1,nz0:nz1))
     endif
     if cont eq 2 or vect eq 2 then begin
       vx3=reform(vx1(kx,ny0:ny1,nz0:nz1))
       vy3=reform(vy1(kx,ny0:ny1,nz0:nz1))
       vz3=reform(vz1(kx,ny0:ny1,nz0:nz1))
     endif
     if cont eq 3 or vect eq 3 then begin
       fx3=reform(ffx1(kx,ny0:ny1,nz0:nz1))
       fy3=reform(ffy1(kx,ny0:ny1,nz0:nz1))
       fz3=reform(ffz1(kx,ny0:ny1,nz0:nz1))
     endif
     if cont eq 4 then b3=reform(alpha1(kx,ny0:ny1,nz0:nz1))
     if cont eq 5 then b3=reform( grad1(kx,ny0:ny1,nz0:nz1))
     if cont eq 6 then b3=reform(  div4(kx,ny0:ny1,nz0:nz1))
     if cont le 3 then begin
       if abso then begin
         case cont of
         0: b3=sqrt(bx3^2+by3^2+bz3^2)
         1: b3=sqrt(cx3^2+cy3^2+cz3^2)
         2: b3=sqrt(vx3^2+vy3^2+vz3^2)
         3: b3=sqrt(fx3^2+fy3^2+fz3^2)
         endcase
         cchr=cchr+'_A'
       endif else begin
         case cont of
         0: b3=bx3
         1: b3=cx3
         2: b3=vx3
         3: b3=fx3
         endcase
         cchr=cchr+'_X'
       endelse
     endif
     case vect of
     0: begin & b1=by3 & b2=bz3 & end
     1: begin & b1=cy3 & b2=cz3 & end
     2: begin & b1=vy3 & b2=vz3 & end
     3: begin & b1=fy3 & b2=fz3 & end
     endcase
     v1=yrange(0)+dely*findgen(nyn)
     v2=zrange(0)+delz*findgen(nzn)
     vchr=vchr+'_YZ'
     xtitle='Y'
     ytitle='Z'
     pchr='X='
     xrange_plt=[yrange(0)-0.5*dely,yrange(1)+0.5*dely]
     yrange_plt=[zrange(0)-0.5*delz,zrange(1)+0.5*delz]
     ratio=float(nzn)/nyn
     end
  3: begin       ; SZ plane
     pos=0.0
     x1h=float(plot1d.x1h)
     y1h=float(plot1d.y1h)
     x2h=float(plot1d.x2h)
     y2h=float(plot1d.y2h)
     szmax=sqrt((x2h-x1h)^2+(y2h-y1h)^2)
     cth=(x2h-x1h)/szmax
     sth=(y2h-y1h)/szmax
     delta=delx
     nsn=fix(szmax/delta)+1
     srange=[0.0,float((nsn-1)*delta)]
     nzn=round((zrange(1)-zrange(0))/delta)+1
     v1=srange(0)+delta*findgen(nsn)
     v2=zrange(0)+delta*findgen(nzn)
     vx=x1h+v1*cth
     vy=y1h+v1*sth
     b1=fltarr(nsn,nzn)
     b2=fltarr(nsn,nzn)
     b3=fltarr(nsn,nzn)
     for j=0,nzn-1 do begin
       for i=0,nsn-1 do begin
         if cont eq 0 or vect eq 0 then  $
            tr_magn,vx(i),vy(i),v2(j),ds3,bx3,by3,bz3,par=0
         if cont eq 1 or vect eq 1 then  $
            tr_magn,vx(i),vy(i),v2(j),ds3,cx3,cy3,cz3,par=1
         if cont eq 2 or vect eq 2 then  $
            tr_magn,vx(i),vy(i),v2(j),ds3,vx3,vy3,vz3,par=2
         if cont eq 3 or vect eq 3 then  $
            tr_magn,vx(i),vy(i),v2(j),ds3,fx3,fy3,fz3,par=3
         if cont ge 4 then begin
            tr_magn,vx(i),vy(i),v2(j),ds3,var,dum,dum,par=cont
            b3(i,j)=var
         endif
         if cont le 3 then begin
           if abso then begin
             case cont of
             0: b3(i,j)=sqrt(bx3^2+by3^2+bz3^2)
             1: b3(i,j)=sqrt(cx3^2+cy3^2+cz3^2)
             2: b3(i,j)=sqrt(vx3^2+vy3^2+vz3^2)
             3: b3(i,j)=sqrt(fx3^2+fy3^2+fz3^2)
             endcase
           endif else begin
             case cont of
             0: b3(i,j)=by3*cth-bx3*sth
             1: b3(i,j)=cy3*cth-cx3*sth
             2: b3(i,j)=vy3*cth-vx3*sth
             3: b3(i,j)=fy3*cth-fx3*sth
             endcase
           endelse
         endif
         case vect of
         0: begin & b1(i,j)=bx3*cth+by3*sth & b2(i,j)=bz3 & end
         1: begin & b1(i,j)=cx3*cth+cy3*sth & b2(i,j)=cz3 & end
         2: begin & b1(i,j)=vx3*cth+vy3*sth & b2(i,j)=vz3 & end
         3: begin & b1(i,j)=fx3*cth+fy3*sth & b2(i,j)=fz3 & end
         endcase
       endfor
     endfor
     if cont le 3 then begin
       if abso then cchr=cchr+'_A' else cchr=cchr+'_N'
     endif
     vchr=vchr+'_SZ'
     xtitle='S'
     ytitle='Z'
     pchr='N='
     xrange_plt=[srange(0)-0.5*delta,srange(1)+0.5*delta]
     yrange_plt=[zrange(0)-0.5*delta,zrange(1)+0.5*delta]
     ratio=float(nzn)/nsn
     end
  endcase
endif else begin                         ; Show horizontal slice
  if n_elements(pos) eq 0 then pos=0.0
  pos=pos>zmin<zmax
  kz=round((pos-zmin)/delz)>0<nz
  pos=zmin+delz*float(kz)
  vchr=vchr(vect)
  cchr=cchr(cont)
  chr=string(kz,format='("MOD (KZ=",i3,")")')
  nx0=round((xrange(0)-xmin)/delx)<(nx-1)  &  xrange(0)=xmin+delx*float(nx0)
  nx1=round((xrange(1)-xmin)/delx)<nx      &  xrange(1)=xmin+delx*float(nx1)
  ny0=round((yrange(0)-ymin)/dely)<(ny-1)  &  yrange(0)=ymin+dely*float(ny0)
  ny1=round((yrange(1)-ymin)/dely)<ny      &  yrange(1)=ymin+dely*float(ny1)
  nxn=nx1-nx0+1
  nyn=ny1-ny0+1
  if cont eq 0 or vect eq 0 then begin
    bx3=bx1(nx0:nx1,ny0:ny1,kz)
    by3=by1(nx0:nx1,ny0:ny1,kz)
    bz3=bz1(nx0:nx1,ny0:ny1,kz)
  endif
  if cont eq 1 or vect eq 1 then begin
    cx3=cx1(nx0:nx1,ny0:ny1,kz)
    cy3=cy1(nx0:nx1,ny0:ny1,kz)
    cz3=cz1(nx0:nx1,ny0:ny1,kz)
  endif
  if cont eq 2 or vect eq 2 then begin
    vx3=vx1(nx0:nx1,ny0:ny1,kz)
    vy3=vy1(nx0:nx1,ny0:ny1,kz)
    vz3=vz1(nx0:nx1,ny0:ny1,kz)
  endif
  if cont eq 3 or vect eq 3 then begin
    fx3=ffx1(nx0:nx1,ny0:ny1,kz)
    fy3=ffy1(nx0:nx1,ny0:ny1,kz)
    fz3=ffz1(nx0:nx1,ny0:ny1,kz)
  endif
  if cont eq 4 then b3=alpha1(nx0:nx1,ny0:ny1,kz)
  if cont eq 5 then b3= grad1(nx0:nx1,ny0:ny1,kz)
  if cont eq 6 then b3=  div4(nx0:nx1,ny0:ny1,kz)
  if cont le 3 then begin
    if abso then begin
      case cont of
      0: b3=sqrt(bx3^2+by3^2+bz3^2)
      1: b3=sqrt(cx3^2+cy3^2+cz3^2)
      2: b3=sqrt(vx3^2+vy3^2+vz3^2)
      3: b3=sqrt(fx3^2+fy3^2+fz3^2)
      endcase
      cchr=cchr+'_A'
    endif else begin
      case cont of
      0: b3=bz3
      1: b3=cz3
      2: b3=vz3
      3: b3=fz3
      endcase
      cchr=cchr+'_Z'
    endelse
  endif
  case vect of
  0: begin & b1=bx3 & b2=by3 & end
  1: begin & b1=cx3 & b2=cy3 & end
  2: begin & b1=vx3 & b2=vy3 & end
  3: begin & b1=fx3 & b2=fy3 & end
  endcase
  v1=xrange(0)+delx*findgen(nxn)
  v2=yrange(0)+dely*findgen(nyn)
  vchr=vchr+'_XY'
  xtitle='X'
  ytitle='Y'
  pchr='Z='
  xrange_plt=[xrange(0)-0.5*delx,xrange(1)+0.5*delx]
  yrange_plt=[yrange(0)-0.5*dely,yrange(1)+0.5*dely]
  ratio=float(nyn)/nxn
endelse
;
;  Limit number of plotted vectors:  *****
;
vv1=v1
vv2=v2
nv1=n_elements(vv1)
nv2=n_elements(vv2)
nvm=max([nv1,nv2])
nvmax=51
if nvm gt nvmax then begin
  fact=(nvm-1)/nvmax+1
  nv1a=nv1/fact
  if fact*nv1a eq nv1 then begin
    k=fact/2
  endif else begin
    k=(nv1-fact*nv1a)/2
    nv1a=nv1a+1
  endelse
  nv1=nv1a
  i1=k+fact*indgen(nv1)
  vv1=vv1(i1)
  b1=b1(i1,*)
  b2=b2(i1,*)
  nv2a=nv2/fact
  if fact*nv2a eq nv2 then begin
    k=fact/2
  endif else begin
    k=(nv2-fact*nv2a)/2
    nv2a=nv2a+1
  endelse
  nv2=nv2a
  i2=k+fact*indgen(nv2)
  vv2=vv2(i2)
  b1=b1(*,i2)
  b2=b2(*,i2)
endif
;
;  Maximum vector length:
;
vmax=sqrt(max(b1^2+b2^2))
if vmax gt 9999. then fmt_vm='f9.1' else fmt_vm='f9.3'
if vmax lt    1. then fmt_vm='f9.5'
if keyword_set(vm) then begin
  print,string(vmax,format='("Actual VMAX=",'+fmt_vm+')')
endif else vm=vmax
vm=vm>1.e-6            ; use specified VM
fact=sqrt(b1^2+b2^2)/vm
length=max(fact)
i=where(fact gt 1.0,count)
if count gt 0 then begin
  b1(i)=b1(i)/fact(i)
  b2(i)=b2(i)/fact(i)
  length=1.0
endif
print,'length=',length
;
;  Maximum contour or greyscale level:
;
cmax=max(abs(b3))
if cmax gt 9999. then fmt_cm='f9.1' else fmt_cm='f9.3'
if cmax lt    1. then fmt_cm='f9.5'
if keyword_set(cm) then begin
  print,string(cmax,format='("Actual CMAX=",'+fmt_cm+')')
endif else cm=cmax
cm=cm>1.e-8           ; use specified CM
if abso then begin
  cmin=min(b3)
  if cmin gt cm then cmin=0.
endif else cmin=-cm
;
;  Position of plot window:
;
if ratio le ratwin then begin
  rat1=1.0
  rat2=ratio/ratwin
endif else begin
  rat1=ratwin/ratio
  rat2=1.0
endelse
;position=[0.52-0.43*rat1,0.52-0.43*rat2,0.52+0.43*rat1,0.52+0.43*rat2]
position=[0.55-0.43*rat1,0.52-0.43*rat2,0.55+0.43*rat1,0.52+0.43*rat2]
;
;  Display greyscale image (greyscale=2):
;
if !d.name eq 'X' then erase
if n_elements(greyscale) eq 0 then greyscale=0
if greyscale eq 2 then begin
  print,string(cmin,cm,  $
        format='("cmin=",'+fmt_cm+',", cmax=",'+fmt_cm+')')
  xoff=fix(!d.x_size* position(0))
  yoff=fix(!d.y_size* position(1))
  xsiz=fix(!d.x_size*(position(2)-position(0)))
  ysiz=fix(!d.y_size*(position(3)-position(1)))
  if !d.name eq 'X' then begin
    tv,grey_scale(congrid(b3,xsiz,ysiz),cmin,cm),xoff,yoff
  endif else begin
    tv,grey_scale(b3,cmin,cm),xoff,yoff,xsiz=xsiz,ysiz=ysiz
  endelse
endif
;
;  Display Quasi-Separatrix Layer (QSL) map (greyscale=3):
;
if greyscale eq 3 then begin
  if ptyp ne 0 then begin
    print,'Error: QSL display requires ptyp=0 (XY).'
  endif else begin
    if n_elements(map_file) eq 0 or n_elements(map_id) eq 0 then begin
      print,'Error: QSL map not available.'
    endif else begin
      file=string(filename,nt,map_id,format='(a,"_",i5.5,"_Q",i1,".sav")')
      if map_file ne dir+file then begin
        print,'Error: QSL map does not belong to current model.'
      endif else begin
        siz=size(map_ndelt)
        map_nx=siz(1)
        map_ny=siz(2)
        i1=round((xrange(0)-map_xmin)/map_int)>0<(map_nx-1)
        i2=round((xrange(1)-map_xmin)/map_int)>0<(map_nx-1)
        j1=round((yrange(0)-map_ymin)/map_int)>0<(map_ny-1)
        j2=round((yrange(1)-map_ymin)/map_int)>0<(map_ny-1)
        if i2 ge i1+2 and j2 ge j1+2 then begin
          map=map_ndelt(i1:i2,j1:j2)
          x1=map_xmin+map_int*(float(i1)-0.5)
          x2=map_xmin+map_int*(float(i2)+0.5)
          y1=map_ymin+map_int*(float(j1)-0.5)
          y2=map_ymin+map_int*(float(j2)+0.5)
          xfact=(position(2)-position(0))/(xrange_plt(1)-xrange_plt(0))
          yfact=(position(3)-position(1))/(yrange_plt(1)-yrange_plt(0))
          posmap=[position(0)+xfact*(x1-xrange_plt(0)),  $
                  position(1)+yfact*(y1-yrange_plt(0)),  $
                  position(0)+xfact*(x2-xrange_plt(0)),  $
                  position(1)+yfact*(y2-yrange_plt(0))]
          xoff=fix(!d.x_size* posmap(0))+1
          yoff=fix(!d.y_size* posmap(1))+1
          xsiz=fix(!d.x_size*(posmap(2)-posmap(0)))
          ysiz=fix(!d.y_size*(posmap(3)-posmap(1)))
          if !d.name eq 'X' then begin
            map=congrid(map,xsiz,ysiz)
            tv,color_scale(map),xoff,yoff,true=1
          endif else begin
            tv,color_scale(map),xoff,yoff,xsiz=xsiz,ysiz=ysiz,true=1
          endelse
        endif else begin
          print,'Warning: QSL map is not in plotted range.'
        endelse
      endelse
    endelse
  endelse
endif
;
;  Plot axes:
;

if keyword_set(epsfile) then begin
plot,xrange_plt,yrange_plt,/nodata,/noerase,     $
  position=position,color=cc1,           $
  xrange=xrange_plt,xstyle=1,xtitle=xtitle,  $
  yrange=yrange_plt,ystyle=1,ytitle=ytitle,xthick=4,ythick=4
endif else begin

plot,xrange_plt,yrange_plt,/nodata,/noerase,     $
  position=position,color=cc1,           $
  xrange=xrange_plt,xstyle=1,xtitle=xtitle,  $
  yrange=yrange_plt,ystyle=1,ytitle=ytitle

endelse
;
;  Plot title:
;
if keyword_set(notitle) then title=' ' else begin
  fmt='(a,"_",i5.5,", ",a,f6.2,", ",a,", VM=",'+fmt_vm   $
                            +',", ",a,", CM=",'+fmt_cm+',", ",a)'
  title=string(filename,nt,pchr,pos,vchr,vm,cchr,cm,chr,format=fmt)
  print,title
  xyouts,0.05,position(3)+0.01,title,/norm,charsize=charsize
endelse
;
;  Vectors:
;
if not keyword_set(novect) then begin
  velovect,b1,b2,vv1,vv2,/noerase,                $
    position=position,length=length,color=cc1,       $
    xrange=xrange_plt,xstyle=1,yrange=yrange_plt,ystyle=1
endif
;
;  Contours (greyscale = 0, 1 or 3):
;
if greyscale le 1 or greyscale eq 3 then begin
  if not keyword_set(nlev) then nlev=5
  fmt='("cmin=",f9.3,", cmax=",f9.3,", clev=",10(f8.3,:,","))'
  if cm gt 9999. then  $
  fmt='("cmin=",f9.1,", cmax=",f9.1,", clev=",10(f8.1,:,","))'
  if cm lt    1. then  $
  fmt='("cmin=",f9.5,", cmax=",f9.5,", clev=",10(f8.5,:,","))'
  if greyscale eq 1 then print,'Contours: solid=positive, dashed=negative' $
                    else print,'Contours: red=positive, green=negative'
  if abso then begin
    clev=cmin+(cm-cmin)*(findgen(nlev)+0.5)/nlev
    print,string(cmin,cm,reverse(clev),format=fmt)
    if greyscale eq 1 then begin   ; B-W contours
      cols=replicate(cc1,nlev)
      lins=replicate(0,nlev)
    endif else begin               ; Color contours
      cols=replicate(red,nlev)
      lins=intarr(nlev)
    endelse
  endif else begin
    clev=cm*2.^(-0.5-findgen(nlev))
    print,string(cmin,cm,clev,format=fmt)
    clev=[-clev,reverse(clev)]
    if greyscale eq 1 then begin   ; B-W contours
      cols=replicate(cc1,2*nlev)
      lins=[replicate(2,nlev),replicate(0,nlev)]
    endif else begin               ; Color contours (red>0,green<0)
    if keyword_set(epsfile) then begin
      cols=[replicate(blue,nlev),replicate(red,nlev)]
      lins=intarr(2*nlev)
    endif else begin
    cols=[replicate(green,nlev),replicate(red,nlev)]
      lins=intarr(2*nlev)
     endelse
    endelse
  endelse
  if keyword_set(epsfile) then begin
  contour,b3,v1,v2,position=position,levels=clev,/noerase,  $
    c_linestyle=lins,c_colors=cols,color=cc1,  $
    xrange=xrange_plt,xstyle=1,yrange=yrange_plt,ystyle=1,c_thick=2
  endif else begin
  contour,b3,v1,v2,position=position,levels=clev,/noerase,  $
    c_linestyle=lins,c_colors=cols,color=cc1,  $
    xrange=xrange_plt,xstyle=1,yrange=yrange_plt,ystyle=1
  endelse
endif
;
;  Connectivity map (greyscale=4):
;
if greyscale eq 4 then begin
  file=string(dir+filename,nt,format='(a,"_",i5.5,"_con.sav")')
  openr,1,file,error=err
  close,1
  if err ne 0 then begin
    print,'Error: Connectivity map not found'
    return
  endif
  restore,file
; print,'blevel=',blevel,'  ','bouter=',bouter
  dx=[-0.5,0.5,0.5,-0.5]
  dy=[-0.5,-0.5,0.5,0.5]
  for n=0,nobipole*2-1 do begin
    for k=0,npix(n)-1 do begin
      if symbol(n,k) eq 4 then  $
        polyfill,xstart(n,k)+dx,ystart(n,k)+dy,color=red,noclip=0
      if symbol(n,k) eq 2 then  $
        polyfill,xstart(n,k)+dx,ystart(n,k)+dy,color=blue,noclip=0
      if symbol(n,k) eq 6 then begin
        if flag(n,k) eq 0 then begin
          polyfill,xstart(n,k)+dx,ystart(n,k)+dy,color=green,noclip=0
        endif else begin
          polyfill,xstart(n,k)+dx,ystart(n,k)+dy,color=yellow,noclip=0
        endelse
      endif
    endfor
  endfor
  for n=0,nobipole*2-1 do begin
    oplot,xboundary(n,0:npts_boundary(n)),  $
          yboundary(n,0:npts_boundary(n)),color=yellow
;   oplot,xboundary2(n,0:npts_boundary2(n)),  $
;         yboundary2(n,0:npts_boundary2(n)),color=yellow,linestyle=1
  endfor
endif
;
;  Show extraction line:
;
if keyword_set(xtra) and n_elements(plot1d) gt 0 then begin
  x1h=float(plot1d.x1h)
  y1h=float(plot1d.y1h)
  x2h=float(plot1d.x2h)
  y2h=float(plot1d.y2h)
  xv =float(plot1d.xv)
  yv =float(plot1d.yv)
  z1v=float(plot1d.z1v)
  z2v=float(plot1d.z2v)
  case ptyp of
  0: oplot,[x1h,x2h],[y1h,y2h],color=yellow,thick=3.
  1: oplot,[xv ,xv ],[z1v,z2v],color=yellow
  2: oplot,[yv ,yv ],[z1v,z2v],color=yellow
  3: begin
       szmax=sqrt((x2h-x1h)^2+(y2h-y1h)^2)
       cth=(x2h-x1h)/szmax
       sth=(y2h-y1h)/szmax
       ss=(xv-x1h)*cth+(yv-y1h)*sth
       oplot,[ss,ss],[z1v,z2v],color=yellow
     end
  endcase
endif
;
;  Show field lines (black/white) and nulls:
;
if n_elements(flid) eq 0 then flid=0
if flid eq 0 then begin
  if n_elements(np) eq 0 then subr_reset
  nnp=np
  if nnp eq 0 then begin
    print,'No starting points'
  endif else begin
    print,'Use starting points in (xp,yp,zp) arrays'
    xxp=xp
    yyp=yp
    zzp=zp
    nnp=np
  endelse
endif else begin
  nnp=0
  flid=fix(flid)>1
  file=string(filename,nt,format='(a,"_",i5.5,"_P")')+  $
       strcompress(flid,/remove_all)
  openr,1,dir+file,error=err
  if err ne 0 then begin
    print,'Error opening file '+file
    close,1
    print,'No starting points'
  endif else begin
    print,'Use starting points in file '+file
    readf,1,nnp
    xxp=fltarr(nnp)
    yyp=fltarr(nnp)
    zzp=fltarr(nnp)
    x0=0.
    y0=0.
    z0=0.
    for n=0,nnp-1 do begin
      readf,1,x0,y0,z0
      xxp(n)=x0
      yyp(n)=y0
      zzp(n)=z0
    endfor
    close,1
  endelse
endelse
;
num=0L
if nnp gt 0 then begin
  if n_elements(flin) eq 0 then flin=-1
  if flin lt 0 then begin
    print,string(nnp,format='("Show all field lines, nnp =",i3)')
    n1=0
    n2=nnp-1
  endif else begin
    flin=flin<(nnp-1)
    print,string(flin,format='("Show one field line, Index =",i3)')
    n1=flin
    n2=flin
  endelse
  if n_elements(smin) eq 0 then smin=-999.
  if n_elements(smax) eq 0 then smax= 999.
  for n=n1,n2 do begin
    trace,xxp(n),yyp(n),zzp(n),smin=smin,smax=smax
    k=n_elements(xtrc)-1
    if num eq 0L then begin
      x2a=xtrc(0:k-1)  &  x2b=xtrc(1:k)
      y2a=ytrc(0:k-1)  &  y2b=ytrc(1:k)
      z2a=ztrc(0:k-1)  &  z2b=ztrc(1:k)
      typ=bytarr(k)
    endif else begin
      x2a=[x2a,xtrc(0:k-1)]  &  x2b=[x2b,xtrc(1:k)]
      y2a=[y2a,ytrc(0:k-1)]  &  y2b=[y2b,ytrc(1:k)]
      z2a=[z2a,ztrc(0:k-1)]  &  z2b=[z2b,ztrc(1:k)]
      typ=[typ,bytarr(k)]
    endelse
    num=num+k
  endfor
  ind=where((x2a-xmin)*(x2b-xmin) gt 0. and (x2a-xmax)*(x2b-xmax) gt 0.)
  kk=where(x2a(ind) lt xmin)
  if kk(0) gt -1 then begin
    index=ind(kk)
    xdel=xmax-xmin
    x2a(index)=x2a(index)+xdel
    x2b(index)=x2b(index)+xdel
  endif
  kk=where(x2a(ind) gt xmax)
  if kk(0) gt -1 then begin
    index=ind(kk)
    xdel=xmax-xmin
    x2a(index)=x2a(index)-xdel
    x2b(index)=x2b(index)-xdel
  endif
  x2a=x2a(ind)  &  x2b=x2b(ind)
  y2a=y2a(ind)  &  y2b=y2b(ind)
  z2a=z2a(ind)  &  z2b=z2b(ind)
  typ=typ(ind)
  num=n_elements(ind)
endif
;
if n_elements(null) eq 0 then null=0
if null gt 0 then begin
  if n_elements(np0) eq 0 then np0=0
  nnp0=np0
  if nnp0 eq 0 then begin
    print,'No nulls'
  endif else begin
    print,'Use nulls in (xp0,yp0,zp0) arrays'
    xxp0=xp0
    yyp0=yp0
    zzp0=zp0
    nnp0=np0
  endelse
endif else nnp0=0
;
if nnp0 gt 0 then begin
  len=2*delx
  dx=[len,0.0,0.0]
  dy=[0.0,len,0.0]
  dz=[0.0,0.0,len]
  for n=0,nnp0-1 do begin
    if xxp0(n) gt xrange(0) and xxp0(n) lt xrange(1) and $
       yyp0(n) gt yrange(0) and yyp0(n) lt yrange(1) and $
       zzp0(n) gt zrange(0) and zzp0(n) lt zrange(1) then begin
      if num eq 0L then begin
        x2a=xxp0(n)-dx  &  x2b=xxp0(n)+dx
        y2a=yyp0(n)-dy  &  y2b=yyp0(n)+dy
        z2a=zzp0(n)-dz  &  z2b=zzp0(n)+dz
        typ=bytarr(3)+1
      endif else begin
        x2a=[x2a,xxp0(n)-dx]  &  x2b=[x2b,xxp0(n)+dx]
        y2a=[y2a,yyp0(n)-dy]  &  y2b=[y2b,yyp0(n)+dy]
        z2a=[z2a,zzp0(n)-dz]  &  z2b=[z2b,zzp0(n)+dz]
        typ=[typ,bytarr(3)+1]
      endelse
      num=num+3
    endif
  endfor
endif
;
if num gt 0 then begin
  case ptyp of
  0: begin
       ind=sort(0.5*(z2a+z2b))
       for nn=0L,num-1 do begin
         n=ind(nn)
         if typ(n) eq 1 then begin
           plots,[x2a(n),x2b(n)],[y2a(n),y2b(n)],/data,  $
            ;     color=cc4,thick=th2bw,noclip=0
            color=cc4,thick=3,noclip=0
         endif else begin
           plots,[x2a(n),x2b(n)],[y2a(n),y2b(n)],/data,  $
                 color=cc2,thick=th2bw,noclip=0
           plots,[x2a(n),x2b(n)],[y2a(n),y2b(n)],/data,  $
                 color=cc1,thick=th1bw,noclip=0
         endelse
       endfor
     end
  1: begin
       ind=sort(-0.5*(y2a+y2b))
       for nn=0L,num-1 do begin
         n=ind(nn)
         if typ(n) eq 1 then begin
           plots,[x2a(n),x2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc4,thick=th2bw,noclip=0
         endif else begin
           plots,[x2a(n),x2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc2,thick=th2bw,noclip=0
           plots,[x2a(n),x2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc1,thick=th1bw,noclip=0
         endelse
       endfor
     end
  2: begin
       ind=sort(0.5*(x2a+x2b))
       for nn=0L,num-1 do begin
         n=ind(nn)
         if typ(n) eq 1 then begin
           plots,[y2a(n),y2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc4,thick=th2bw,noclip=0
         endif else begin
           plots,[y2a(n),y2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc2,thick=th2bw,noclip=0
           plots,[y2a(n),y2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc1,thick=th1bw,noclip=0
         endelse
       endfor
     end
  3: begin
       tmp= (x2a-x1h)*cth+(y2a-y1h)*sth
       y2a=-(x2a-x1h)*sth+(y2a-y1h)*cth
       x2a=tmp
       tmp= (x2b-x1h)*cth+(y2b-y1h)*sth
       y2b=-(x2b-x1h)*sth+(y2b-y1h)*cth
       x2b=tmp
       ind=sort(-0.5*(y2a+y2b))
       for nn=0L,num-1 do begin
         n=ind(nn)
         if typ(n) eq 1 then begin
           plots,[x2a(n),x2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc4,thick=th2bw,noclip=0
         endif else begin
           plots,[x2a(n),x2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc2,thick=th2bw,noclip=0
           plots,[x2a(n),x2b(n)],[z2a(n),z2b(n)],/data,  $
                 color=cc1,thick=th1bw,noclip=0
         endelse
       endfor
     end
  endcase
endif
;
;  Close Postscript file:
;
if keyword_set(epsfile) then begin
  device,/close
  print,'Plot saved on '+epsfile
  set_plot,'X'
  hexa_colour
endif
;
print,'-----'
end

;------------------------------------------------------------------
;  Create and Realize 3D Display Window:
;------------------------------------------------------------------
pro PROG_PLOT3D
@hexa.blk
@hexa.trc
@hexa.col
common wid,wid_main
;
;  Setup W3D structure:
;
w3d={w3d, draw: 0L,  alt: 0L, azim: 0L, wdth: 0L, xcen: 0L, ycen: 0L, $
          zcen: 0L, grey: 0L, cmax: 0L, nlev: 0L, bbox: 0L, flid: 0L, $
          flin: 0L, smin: 0L, smax: 0L, flbw: 0L, snip: 0L, bars: 0L, $
          snum: 0L, zbse: 0L, zstr: 0L, sdel: 0L, nwin: 0L, $
          alt2: 0L, azm2: 0L, wdt2: 0L,  nt1: 0L,  nt2: 0L}


;           w3d={w3d, draw: 0L,  alt: 0L, azim: 0L, wdth: 0L, xcen: 0L, ycen: 0L, $
;           zcen: 0L, grey: 0L, cmax: 0L, nlev: 0L, bbox: 0L, flid: 0L, $
;           flin: 0L, smin: 0L, smax: 0L, flbw: 0L, snip: 0L, bars: 0L, $
;           snum: 0L, zbse: 0L, zstr: 0L, sdel: 0L, null: 0L, nwin: 0L, $
;           alt2: 0L, azm2: 0L, wdt2: 0L,  nt1: 0L,  nt2: 0L}
;
;  Create base widget:
;
base=widget_base(title='3D Display',group_leader=wid_main(0),  $
                 xoffset=400,yoffset=0,/column,uvalue=w3d)
;
;  First row contains 2D image display:
;
row1=widget_base(base,/row)
w3d.draw=widget_draw(row1,/button_events,uvalue='DRAW',   $
                           x_scroll_size=800,xsize=820,   $
                           y_scroll_size=700,ysize=720)
;
;  Second row contains image control widgets:
;
row2=widget_base(base,/row)
;
;  1) First group contains viewing parameters, VIEW and QUIT buttons:
;
acol=widget_base(row2,/column)
grp=widget_base(acol,/row)
alt=30.
azim=30.
width=1.5*xmax
w3d.alt =cw_field(grp,title='ALT:' ,value=alt ,/floating,xsize=5)
w3d.azim=cw_field(grp,title='AZIM:',value=azim,/floating,xsize=5)
grp=widget_base(acol,/row)
w3d.wdth=cw_field(grp,title='Width:',value=width,/floating,xsize=5)
item=widget_button(grp,value='Rs',uvalue='RESET-LLW')
item=widget_button(acol,value='VIEW',uvalue='VIEW')
item=widget_button(acol,value='F-L Reset (ID=0)',uvalue='RESET-FLN')
item=widget_button(acol,value='View in Window 0',uvalue='VIEW0')
item=widget_button(acol,value='EPS File',uvalue='EPSF')
item=widget_button(acol,value='QUIT',uvalue='QUIT')
;
;  2) Second group for focus position:
;
acol=widget_base(row2,/column)
item=widget_label(acol,value='Focus:')
xpos=0.5*(xmin+xmax)
ypos=0.5*(ymin+ymax)
w3d.xcen=cw_field(acol,title='X:',value=xpos,/floating,xsize=5)
w3d.ycen=cw_field(acol,title='Y:',value=ypos,/floating,xsize=5)
w3d.zcen=cw_field(acol,title='Z:',value=zmin,/floating,xsize=5)
item=widget_button(acol,value='Reset',uvalue='RESET-FOCUS')
;
;  3) Third group controls photospheric magnetogram, connectivity map,
;     background image, field lines, bars (for F-L dips or fibrils),
;     bounding box, F-L selection params, coronal loop selection params:
;
acol=widget_base(row2,/column)
grp=widget_base(acol,/row)
w3d.grey=cw_bgroup(grp,['CC','BW','Grey','No'],column=4,  $
                  /exclusive,/no_release,label_left='Surface:',  $
                  set_value=0,uvalue='GREYSCALE')
w3d.cmax=cw_field(grp,title='Max:',value=0.,/floating,xsize=8)
w3d.nlev=cw_field(grp,title='Nlev:',value=5,/integer,xsize=3)
grp=widget_base(acol,/row)
w3d.flid=cw_field(grp,title='F-L ID:',value=0,/integer,xsize=3)
w3d.flin=cw_field(grp,title='Index:',value=-1,/integer,xsize=3)
w3d.smin=cw_field(grp,title='SMin:',value=-999.,/float,xsize=6)
w3d.smax=cw_field(grp,title='SMax:',value= 999.,/float,xsize=6)
w3d.flbw=cw_bgroup(grp,['Col','BW'],column=2,  $
                  /exclusive,/no_release,label_left=' ', $
                  set_value=0,uvalue='FLBW')
;w3d.snip=cw_bgroup(grp,['No','Yes'],column=2,  $
;                  /exclusive,/no_release,label_left='Snip:', $
;                  set_value=0,uvalue='SNIP')
grp=widget_base(acol,/row)
w3d.snum=cw_field(grp,title='F-L Select Num:',value=1,  $
                  /integer,xsize=3)
w3d.zbse=cw_field(grp,title='zbase:' ,value=zmin,/floating,xsize=5)
w3d.zstr=cw_field(grp,title='zstart:',value=zmin,/floating,xsize=5)
w3d.sdel=cw_field(grp,title='Del:'   ,value=sdel,/floating,xsize=5)
grp=widget_base(acol,/row)
;w3d.bars=cw_bgroup(grp,['No','Dips','Fibrils'],column=3,  $
w3d.bars=cw_bgroup(grp,['No','Dips'],column=2,  $
                  /exclusive,/no_release,label_left='Bars:', $
                  set_value=0,uvalue='BARS')
;w3d.null=cw_bgroup(grp,['No','Yes'],column=2,  $
 ;                  /exclusive,/no_release,label_left='Nulls:',  $
 ;                  set_value=0,uvalue='NULL')
;item=widget_button(grp,value='Fan/Spine',uvalue='NSEL')
w3d.bbox=cw_bgroup(grp,['No','Yes'],column=2,  $
                  /exclusive,/no_release,label_left='Box:', $
                  set_value=1,uvalue='BOX')
grp=widget_base(acol,/row)
alt2=alt
azim2=azim+90.
wdth2=width
item=widget_button(grp,value='Rotate',uvalue='ROT')
w3d.alt2=cw_field(grp,title='ALT2:' ,value=alt2 ,/floating,xsize=5)
w3d.azm2=cw_field(grp,title='AZIM2:',value=azim2,/floating,xsize=5)
w3d.wdt2=cw_field(grp,title='WDTH2:',value=wdth2,/floating,xsize=5)
item=widget_button(grp,value='TimeSeq',uvalue='TIM')
nt1=nrelax
nt2=nrelax+nmajor
w3d.nt1=cw_field(grp,title='nt1:',value=nt1,/integer,xsize=4)
w3d.nt2=cw_field(grp,title='nt2:',value=nt2,/integer,xsize=4)
;
;  Realize widgets and register with XMANAGER:
;
widget_control,/realize,base
widget_control,w3d.draw,get_value=nwin  ; get IDL window index
w3d.nwin=nwin
widget_control,base,set_uvalue=w3d      ; store W3D in UVALUE of base widget
plot3d,nwin=w3d.nwin,alt=alt,azim=azim,width=width,/box
xmanager,'PLOT3D',base,/no_block    ; /just_reg
end

;-------------------------------------------------------------------------
;  Process Events from 3D Display Window:
;-------------------------------------------------------------------------
pro PLOT3D_EVENT,ev
@hexa.blk
@hexa.trc
@hexa.col
;
;  Get UVALUE of event:
;
widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=w3d
;
if uval eq 'QUIT' then begin
  widget_control,ev.top,/destroy
  return
endif
;
;  Reset viewing window (ALT, AZIM and WIDTH):
;
if uval eq 'RESET-LLW' then begin
  widget_control,w3d.alt ,set_value=30.0
  widget_control,w3d.azim,set_value=30.0
  widget_control,w3d.wdth,set_value=1.5*xmax
  return
endif
;
;  Reset focus position:
;
if uval eq 'RESET-FOCUS' then begin
  xpos=0.5*(xmin+xmax)
  ypos=0.5*(ymin+ymax)
  widget_control,w3d.xcen,set_value=xpos
  widget_control,w3d.ycen,set_value=ypos
  widget_control,w3d.zcen,set_value=zmin
  return
endif
widget_control,w3d.zbse,set_value=zbase
widget_control,w3d.zstr,set_value=zstart
;
;  Select nearest null:
;
if uval eq 'NSEL' then begin
  if np0 eq 0 then begin
    POP_UP,'Error','No nulls selected'
    return
  endif
  if w3d.nwin ne !d.window then begin
    POP_UP,'Error','First VIEW this window before selecting null'
    return
  endif
  widget_control,w3d.alt ,get_value=alt
  widget_control,w3d.azim,get_value=azim
  subr_rotate,alt,azim,xp0,yp0,zp0,xv0,yv0,zv0
  print,'Click LEFT in PLOT3D to select a NULL point (black cross)'
  cursor,xv1,yv1,3,/data
  if !err eq 4 then return
  dstmin=min(sqrt((xv0-xv1)^2+(yv0-yv1)^2),nn)
  print,string(nn,xp0(nn),yp0(nn),zp0(nn),  $
        format='("nn=",i5,",  xp0=",f7.3,",  yp0=",f7.3,",  zp0=",f7.3)')
  prog_null_launch,nn
  uval='VIEW'
endif
;
;  Show 3D Display:
;
if uval eq 'VIEW' or uval eq 'VIEW0' or uval eq 'EPSF' or  $
   uval eq 'RESET-FLN' or uval eq 'RESET-CLP' then begin
  if uval eq 'RESET-FLN' then begin
    subr_reset
;    widget_control,w3d.bars,set_value=0
  endif
  if uval eq 'EPSF' then epsfile=dir+'idl.eps' else epsfile=''
  if uval eq 'VIEW0' then nwin=0 else nwin=w3d.nwin
  widget_control,w3d.alt ,get_value=alt
  widget_control,w3d.azim,get_value=azim
  widget_control,w3d.wdth,get_value=width
  widget_control,w3d.xcen,get_value=xcen
  widget_control,w3d.ycen,get_value=ycen
  widget_control,w3d.zcen,get_value=zcen
  widget_control,w3d.grey,get_value=greyscale
  widget_control,w3d.cmax,get_value=cm
  widget_control,w3d.nlev,get_value=nlev
  widget_control,w3d.bbox,get_value=box
  widget_control,w3d.flid,get_value=flid
  widget_control,w3d.flin,get_value=flin
  widget_control,w3d.smin,get_value=smin
  widget_control,w3d.smax,get_value=smax
  widget_control,w3d.flbw,get_value=flbw
;  widget_control,w3d.snip,get_value=snip
  widget_control,w3d.zbse,get_value=zbase
  widget_control,w3d.bars,get_value=bars
 ; widget_control,w3d.null,get_value=null
  plot3d,alt=alt,azim=azim,width=width,xcen=xcen,ycen=ycen,zcen=zcen,  $
         greyscale=greyscale,cm=cm,nlev=nlev,box=box,  $
         flid=flid,flin=flin,smin=smin,smax=smax,flbw=flbw,zbase=zbase,  $
         bars=bars,null=null,nwin=nwin,epsfile=epsfile
  return
endif
;
;  Rotate 3D model:
;
if uval eq 'ROT' then begin
  nwin=w3d.nwin
  widget_control,w3d.alt ,get_value=alt1
  widget_control,w3d.azim,get_value=azim1
  widget_control,w3d.wdth,get_value=wdth1
  widget_control,w3d.xcen,get_value=xcen
  widget_control,w3d.ycen,get_value=ycen
  widget_control,w3d.zcen,get_value=zcen
  widget_control,w3d.grey,get_value=greyscale
  widget_control,w3d.cmax,get_value=cm
  widget_control,w3d.nlev,get_value=nlev
  widget_control,w3d.bbox,get_value=box
  widget_control,w3d.flid,get_value=flid
  widget_control,w3d.flin,get_value=flin
  widget_control,w3d.smin,get_value=smin
  widget_control,w3d.smax,get_value=smax
  widget_control,w3d.flbw,get_value=flbw
  widget_control,w3d.zbse,get_value=zbase
  widget_control,w3d.bars,get_value=bars
 ; widget_control,w3d.null,get_value=null
  widget_control,w3d.alt2,get_value=alt2
  widget_control,w3d.azm2,get_value=azim2
  widget_control,w3d.wdt2,get_value=wdth2
  alnw1=alog(wdth1)
  alnw2=alog(wdth2)
  dalt_max =0.5   ; step size in degrees
  dazim_max=0.5
  dlnw_max=0.01
  nframe=max([round((alt2 -alt1 )/dalt_max )+1,  $
              round((azim2-azim1)/dazim_max)+1,  $
              round((alnw2-alnw1)/dlnw_max )+1,3])
  dalt =(alt2 -alt1 )/(nframe-1)
  dazim=(azim2-azim1)/(nframe-1)
  dlnw =(alnw2-alnw1)/(nframe-1)
  answer=''
  ;read,prompt='Close XINTERANIMATE widget? (y/n) ',answer

  ;if answer eq 'y' then xinteranimate,/close
  xinteranimate,set=[820,520,nframe],/showload,/track,/cycle, $
                title='Plot3d Movie'
  for n=0,nframe-1 do begin
    alt =alt1 +dalt *n
    azim=azim1+dazim*n
    alnw=alnw1+dlnw *n
    width=exp(alnw)
    plot3d,alt=alt,azim=azim,width=width,xcen=xcen,ycen=ycen,zcen=zcen,  $
           greyscale=greyscale,cm=cm,nlev=nlev,box=box,flid=flid,flin=flin, $
           smin=smin,smax=smax,flbw=flbw,zbase=zbase,bars=bars,null=null,  $
           nwin=nwin
    xinteranimate,frame=n,window=[nwin,0,100,820,520]
  endfor
  xinteranimate,20.
endif
;
;  Time sequence:
;
if uval eq 'TIM' then begin
  nwin=w3d.nwin
  widget_control,w3d.alt ,get_value=alt
  widget_control,w3d.azim,get_value=azim
  widget_control,w3d.wdth,get_value=width
  widget_control,w3d.xcen,get_value=xcen
  widget_control,w3d.ycen,get_value=ycen
  widget_control,w3d.zcen,get_value=zcen
  widget_control,w3d.grey,get_value=greyscale
  widget_control,w3d.cmax,get_value=cm
  widget_control,w3d.nlev,get_value=nlev
  widget_control,w3d.bbox,get_value=box
  widget_control,w3d.flid,get_value=flid
  widget_control,w3d.flin,get_value=flin
  widget_control,w3d.smin,get_value=smin
  widget_control,w3d.smax,get_value=smax
  widget_control,w3d.flbw,get_value=flbw
  widget_control,w3d.zbse,get_value=zbase
  widget_control,w3d.bars,get_value=bars
 ; widget_control,w3d.null,get_value=null
  widget_control,w3d.nt1 ,get_value=nt1
  widget_control,w3d.nt2 ,get_value=nt2
  nframe=nt2-nt1+1
  answer=''
  ;read,prompt='Close XINTERANIMATE widget? (y/n) ',answer
  ;if answer eq 'y' then xinteranimate,/close
  xinteranimate,set=[820,520,5*nframe],/showload,/track, $
                title='Plot3d Movie'
  for n=nt1,nt2 do begin
    nt=long(n)
    fullname=string(filename,nt,format='(a,"_",i5.5)')+extension
    prog_restore,fullname
    plot3d,alt=alt,azim=azim,width=width,xcen=xcen,ycen=ycen,zcen=zcen,  $
           greyscale=greyscale,cm=cm,nlev=nlev,box=box,flid=flid,flin=flin, $
           smin=smin,smax=smax,flbw=flbw,zbase=zbase,bars=bars,null=null,  $
           nwin=nwin
    xinteranimate,frame=5*(nt-nt1),window=[nwin,0,100,820,520]
  endfor
  xinteranimate,20.
endif
;
;  Trace field lines (click in DRAW widget):
;
if uval eq 'DRAW' then begin
  if ev.press ne 1 then return   ; only accept Left Button events
  if w3d.nwin ne !d.window then begin
    POP_UP,'Error','First VIEW this window before selecting field line(s)'
    return
  endif
  widget_control,w3d.alt ,get_value=alt      ; viewing angles
  widget_control,w3d.azim,get_value=azim
  widget_control,w3d.smin,get_value=smin
  widget_control,w3d.smax,get_value=smax
  widget_control,w3d.zstr,get_value=zz0      ; starting height
  widget_control,w3d.snum,get_value=snum     ; number of points along LOS
  snum=snum>1<99
  widget_control,w3d.sdel,get_value=sdel      ; spacing parameter
  sdel=sdel>0.01
  widget_control,w3d.sdel,set_value=sdel      ; spacing parameter
  print,string(alt,azim,  $
        format='("alt=",f7.3," azim=",f8.3," (deg)")')
  xv=(float(ev.x)/!d.x_size-!x.s(0))/!x.s(1)
  yv=(float(ev.y)/!d.y_size-!y.s(0))/!y.s(1)
  subr_invert1,alt,azim,xv,yv,xx0,yy0,ind,zlevel=zz0
  if ind(0) eq -1 then begin
    print,'Position not in box:'
    print,string(xx0,yy0,zz0,  $
          format='("xx0=",f7.3,"  yy0=",f7.3,"  zz0=",f7.3)')
    return
  endif
  theta=(90.0-alt)/!radeg
  azimr=azim/!radeg
  if n_elements(np) eq 0 then subr_reset
  for n=0,snum-1 do begin
    pos=sdel*n
    xx=xx0-pos*sin(theta)*sin(azimr)
    yy=yy0-pos*sin(theta)*cos(azimr)
    zz=zz0+pos*cos(theta)
    print,string(n,xx,yy,zz,  $
          format='("n=",i2,"  xx=",f7.3,"  yy=",f7.3,"  zz=",f7.3)')
    if np eq 0 then begin
      xp(0)=xx
      yp(0)=yy
      zp(0)=zz
    endif else begin
      xp=[xp,xx]
      yp=[yp,yy]
      zp=[zp,zz]
    endelse
    np=np+1
    trace,xx,yy,zz,zbase=zbase,smin=smin,smax=smax
    num=n_elements(xtrc)
    subr_rotate,alt,azim,xtrc,ytrc,ztrc,xs2,ys2,zs2
    oplot,xs2,ys2,color=black
  endfor
endif
end

;------------------------------------------------------------------
;  Plot magnetic field lines in 3D:
;------------------------------------------------------------------
pro plot3d,alt=alt,azim=azim,width=width,xcen=xcen,ycen=ycen,zcen=zcen, $
           greyscale=greyscale,cm=cm,nlev=nlev,flid=flid,flin=flin,  $
           smin=smin,smax=smax,flbw=flbw,zbase=zbase,bars=bars, $
           null=null,box=box,nwin=nwin,epsfile=epsfile
@hexa.blk
@hexa.trc
@hexa.col
;
;  Setup plot device:
;
if keyword_set(epsfile) then begin
  set_plot,'PS'
  nxv=1200                  ; number of pixels
  nyv=1200
  xsiz=0.015*nxv           ; image size in cm
  ysiz=0.015*nyv
  device,/encapsulated,/portrait,/color,bits=8,  $
         filename=epsfile,xsize=xsiz,ysize=ysiz
  hexa_colour
  cc2=white  ; whiteout (background) color
  cc1=black  ; drawing color
  cc4=red    ; nulls
  th2bx=9.   ; whiteout line thickness for box
  th1bx=2.   ; drawing line thickness for box and contours
  th2bw=12.  ; whiteout line thickness for field lines (black-white mode)
  th1bw=5.   ; drawing line thickness for field lines (black-white mode)
  thcol=7.   ; drawing line thickness for field lines (color mode)
endif else begin
  set_plot,'X'
  if n_elements(nwin) eq 0 then nwin=0
  if nwin gt 0 then wset,nwin else window,0,xsize=820,ysize=720
  nxv=!d.x_size
  nyv=!d.y_size
  hexa_colour
  cc2=white  ; whiteout (background) color
  cc1=black  ; drawing color
;  cc2=black  ; whiteout color
;  cc1=white  ; drawing color
  cc4=black  ; nulls
  th2bx=5.   ; whiteout line thickness for box
  th1bx=1.   ; drawing line thickness for box and contours
  th2bw=6.   ; whiteout line thickness for field lines (black-white mode)
  th1bw=2.   ; drawing line thickness for field lines (black-white mode)
  thcol=3.   ; drawing line thickness for field lines (color mode)
endelse
erase,cc2
;
;  Set up ALT-AZIM viewing angles:
;
if n_elements(alt) eq 0 then alt=0.0
alt=float(alt)>0.<90.
if n_elements(azim) eq 0 then azim=0.0
azim=float(azim)>(-180.)<180.
print,string(alt,azim,format='("Observer: alt=",f8.3,", azim=",f8.3," deg")')
;
;  Determine focal point (put at center of viewing window):
;
if n_elements(xcen) eq 0 then xcen=0.5*(xmin+xmax)
if n_elements(ycen) eq 0 then ycen=0.5*(ymin+ymax)
if n_elements(zcen) eq 0 then zcen=zmin
print,string(xcen,ycen,zcen,  $
      format='("xcen=",f6.1," ycen=",f6.1," zcen=",f6.1)')
subr_rotate,alt,azim,xcen,ycen,zcen,xvcen,yvcen,zvcen
;
;  Determine viewing window in solar coordinates:
;
if n_elements(width) eq 0 then width=1.5*xmax
width=float(width)
ratwin=float(!d.y_size)/!d.x_size
xvmin=xvcen-0.5*width
xvmax=xvcen+0.5*width
yvmin=yvcen-0.5*width*ratwin
yvmax=yvcen+0.5*width*ratwin
pix=width/nxv             ; pixel size
;
;  Display magnetogram as a greyscale image (greyscale=2):
;
if n_elements(cm) eq 0 then cm=0.
if cm le 0. then cm=max(abs(bz1(*,*,0)))
if n_elements(greyscale) eq 0 then greyscale=0
if greyscale eq 2 then begin
  plot,[0,1],[0,1],/nodata,xstyle=5,xrange=[xvmin,xvmax],  $
                           ystyle=5,yrange=[yvmin,yvmax],  $
                           position=[0.,0.,1.,1.]
  xv=(xvmin+pix*findgen(nxv))#(fltarr(nyv)+1.0)
  yv=(fltarr(nxv)+1.0)#(yvmin+pix*findgen(nyv))
  subr_invert1,alt,azim,xv,yv,x,y,ind
  if ind(0) gt -1 then begin
    aa=fltarr(nxv,nyv)
    if cc2 eq black then aa=aa-cm
    if cc2 eq white then aa=aa+cm
    aa(ind)=interpolate(bz1(*,*,0),(x(ind)-xmin)/delx-0.5,  $
                                   (y(ind)-ymin)/dely-0.5)
    img=grey_scale(aa,-cm,cm)
    if !d.name eq 'X' then tv,img else tv,img,xsiz=1000*xsiz,ysiz=1000*ysiz
  endif else print,'Magnetogram not visible'
endif
;
;  Draw magnetogram contours (greyscale=0 or 1):
;
if greyscale le 1 then begin
  if n_elements(nlev) eq 0 then nlev=5
  blev=cm*2.^(-0.5-findgen(nlev))
  print,string(cm,blev,format='("bmax=",f8.3,", blev=",10(f8.3,:,","))')
  blev=[-blev,reverse(blev)]
  if greyscale eq 0 then begin    ; Color contours
    cols=[replicate(green,nlev),replicate(red,nlev)]  ; red>0, green<0
    lins=intarr(2*nlev)
  endif else begin                ; B-W contours
    cols=replicate(cc1,2*nlev)
    lins=[replicate(2,nlev),replicate(0,nlev)]        ; solid>0, dashed<0
  endelse
  xa=xmin+delx*findgen(nx+1)
  ya=ymin+dely*findgen(ny+1)
  contour,bz1(*,*,0),xa,ya,levels=blev,  $
          path_info=info,path_xy=xy,/path_data_coords
  plot,[0,1],[0,1],/nodata,/noerase,  $
                   xstyle=5,xrange=[xvmin,xvmax],  $
                   ystyle=5,yrange=[yvmin,yvmax],  $
                   position=[0.,0.,1.,1.]
  for ii=0,n_elements(info)-1 do begin
    lev=info(ii).level
    num=info(ii).n
    ind=info(ii).offset+[indgen(num),0]
    x=reform(xy(0,ind))
    y=reform(xy(1,ind))
    z=fltarr(num+1)
    subr_rotate,alt,azim,x,y,z,xv,yv,zv
;    oplot,xv,yv,color=cols(lev),linestyle=lins(lev),thick=th1bx
    vis=bytarr(num+1)+1
    edge=where(x lt (xmin+0.3*delx) or x gt (xmax-0.3*delx) or  $
               y lt (ymin+0.3*dely) or y gt (ymax-0.3*dely),count)
    if count gt 0 then vis(edge)=0
    ivis=where(vis(0:num-1) and vis(1:num),count)
    if count gt 0 then begin
      for n=0,count-1 do begin
        i=ivis(n)
        plots,[xv(i),xv(i+1)],[yv(i),yv(i+1)],/data, $  ; noclip=0,  $
              color=cols(lev),linestyle=lins(lev),thick=th1bx
      endfor
    endif
  endfor
endif
;
;  Draw border of magnetogram:
;
if keyword_set(box) then begin
  x=[xmin,xmax,xmax,xmin,xmin]
  y=[ymin,ymin,ymax,ymax,ymin]
  z=[zmin,zmin,zmin,zmin,zmin]
  subr_rotate,alt,azim,x,y,z,xv,yv,zv
  oplot,xv,yv,color=cc1
endif
;
;  Define starting points for field lines:
;
if n_elements(flid) eq 0 then flid=0
if flid eq 0 then begin
  if n_elements(np) eq 0 then subr_reset
  nnp=np
  if nnp eq 0 then begin
    print,'No starting points'
  endif else begin
    print,'Use starting points in (xp,yp,zp) arrays'
    xxp=xp
    yyp=yp
    zzp=zp
    nnp=np
  endelse
endif else begin
  nnp=0
  flid=fix(flid)>1
  file=string(filename,nt,format='(a,"_",i5.5,"_P")')+  $
       strcompress(flid,/remove_all)
  openr,1,dir+file,error=err
  if err ne 0 then begin
    print,'Error opening file '+file
    close,1
    print,'No starting points'
  endif else begin
    print,'Use starting points in file '+file
    readf,1,nnp
    xxp=fltarr(nnp)
    yyp=fltarr(nnp)
    zzp=fltarr(nnp)
    x0=0.
    y0=0.
    z0=0.
    for n=0,nnp-1 do begin
      readf,1,x0,y0,z0
      xxp(n)=x0
      yyp(n)=y0
      zzp(n)=z0
    endfor
    close,1
  endelse
endelse
;
;  Assign colors to starting points:
;
if nnp gt 0 then begin
  ncolors=n_elements(colors)
  ccp=colors(indgen(nnp) mod ncolors)     ; cycle through all colors
endif
;
;  Trace field lines in 3D magnetic model:
;
if nnp gt 0 then begin
  if n_elements(flin) eq 0 then flin=-1
  flin=flin>(-1)<(nnp-1)
  if flin eq -1 then begin
    print,string(nnp,format='("Show all field lines (nnp =",i3,")")')
    n1=0   &  n2=nnp-1
  endif else begin
    print,string(flin,format='("Show selected field line, Index =",i3)')
    n1=flin  &  n2=flin
  endelse
  if n_elements(zbse) eq 0 then zbse=zmin
  if n_elements(smin) eq 0 then smin=-999.
  if n_elements(smax) eq 0 then smax= 999.
;
;  Trace field lines and break them up into sections:
;
  for n=n1,n2 do begin
    trace,xxp(n),yyp(n),zzp(n),zbase=zbase,smin=smin,smax=smax
    k=n_elements(xtrc)-1
    if n eq n1 then begin
      x2a=xtrc(0:k-1)  &  x2b=xtrc(1:k)
      y2a=ytrc(0:k-1)  &  y2b=ytrc(1:k)
      z2a=ztrc(0:k-1)  &  z2b=ztrc(1:k)
      ccc=replicate(ccp(n),k)
    endif else begin
      x2a=[x2a,xtrc(0:k-1)]  &  x2b=[x2b,xtrc(1:k)]
      y2a=[y2a,ytrc(0:k-1)]  &  y2b=[y2b,ytrc(1:k)]
      z2a=[z2a,ztrc(0:k-1)]  &  z2b=[z2b,ztrc(1:k)]
      ccc=[ccc,replicate(ccp(n),k)]
    endelse
  endfor
  num2=n_elements(x2a)
  num1=num2
;
;  If keyword SNIP is set, then force all field lines to lie within
;  the computational domain (i.e., remove any sections that cross the
;  side boundaries, and modify x and y coordinates):
;
  if keyword_set(snip)  then begin
    ind=where((x2a-xmin)*(x2b-xmin) gt 0. and (x2a-xmax)*(x2b-xmax) gt 0.)
    num1=n_elements(ind)
    xdel=xmax-xmin
    kk=where(x2a(ind) lt xmin)
    if kk(0) gt -1 then begin
      index=ind(kk)
      x2a(index)=x2a(index)+xdel
      x2b(index)=x2b(index)+xdel
    endif
    kk=where(x2a(ind) gt xmax)
    if kk(0) gt -1 then begin
      index=ind(kk)
      x2a(index)=x2a(index)-xdel
      x2b(index)=x2b(index)-xdel
    endif
  endif else ind=lindgen(num1)
endif else begin
  num1=0L
  num2=0L
endelse
;
;  Plot nulls as 3D crosses:
;
if n_elements(null) eq 0 then null=0
if null gt 0 then begin
  if n_elements(np0) eq 0 then np0=0
  nnp0=np0
  if nnp0 eq 0 then begin
    print,'No nulls'
  endif else begin
    print,'Use nulls in (xp0,yp0,zp0) arrays'
    xxp0=xp0
    yyp0=yp0
    zzp0=zp0
    nnp0=np0
  endelse
endif else nnp0=0
;
if nnp0 gt 0 then begin
  len=2*delx
  dx=[len,0.0,0.0]
  dy=[0.0,len,0.0]
  dz=[0.0,0.0,len]
  for n=0,nnp0-1 do begin
    if num2 eq 0L then begin
      x2a=xxp0(n)-dx  &  x2b=xxp0(n)+dx
      y2a=yyp0(n)-dy  &  y2b=yyp0(n)+dy
      z2a=zzp0(n)-dz  &  z2b=zzp0(n)+dz
      ccc=replicate(cc4,3)
      ind=lindgen(3)
    endif else begin
      x2a=[x2a,xxp0(n)-dx]  &  x2b=[x2b,xxp0(n)+dx]
      y2a=[y2a,yyp0(n)-dy]  &  y2b=[y2b,yyp0(n)+dy]
      z2a=[z2a,zzp0(n)-dz]  &  z2b=[z2b,zzp0(n)+dz]
      ccc=[ccc,replicate(cc4,3)]
      ind=[ind,lindgen(3)+num2]
    endelse
    num1=num1+3
    num2=num2+3
  endfor
endif
;
;  Define bars for dips/fibrils near filament path:
;
if n_elements(bars) eq 0 then bars=0
if bars gt 0 then begin
;
;  a) Field-line dips:
;
  if bars eq 1 then begin
    i1=1   &  i2=nx-1       ; range in x for finding dips (cells)
    j1=10  &  j2=ny-10      ; range in y for finding dips (cells)
    nzmin=1                 ; starting height to search for dips
    nzmax=60<nz             ; maximum search height,
    amin=0.  ; -0.01        ; minimum angle at +/-1 pixel from dip
                            ; (small negative for slight arch, was -0.02)
    print,'Find field-line dips'
    print,string(nzmin,nzmax,amin,  $
          format='("nzmin=",i2,", nzmax=",i2,", amin=",f5.2)')
    print,'        LtBlueFr          Nx    Ny    Nz'
    ndips_tot=0L
    nmx=nzmax
    for n=nzmax,nzmin,-1 do begin    ; loop over heights
      fblue=float(nmx-n)/(nmx-nzmin)
;        light-to-dark blue for dips:
      if !d.name eq 'X' then cc3=256L*(round(255.*fblue)+256L*255)  $
                        else cc3=ncol+round((nfbl-1)*fblue)
      ndips=intarr(3)
      if n eq 0     then m1=2 else m1=0
      if n eq nzmax then m2=1 else m2=2
      for m=m1,m2 do begin
        case m of
        0: begin  ; x-ribs
           index=where(bz1(i1:i2-1,j1:j2,n)*bz1(i1+1:i2,j1:j2,n) lt 0.0,nzero)
           if nzero gt 0 then begin
             x0=fltarr(nzero)  &  y0=fltarr(nzero)  &  z0=fltarr(nzero)
             dx=fltarr(nzero)  &  dy=fltarr(nzero)
             bm0=fltarr(nzero)
             dbz=fltarr(nzero)
             for kz=0,nzero-1 do begin
               i=i1+(index(kz) mod (i2-i1))
               j=j1+index(kz)/(i2-i1)
               fx0=bz1(i,j,n)/(bz1(i,j,n)-bz1(i+1,j,n))
               fx1=1.0-fx0
               x0(kz)=xmin+delx*(float(i)+fx0)
               y0(kz)=ymin+dely*float(j)
               z0(kz)=zmin+delz*float(n)
               bxdip=fx1*bx1(i,j,n)+fx0*bx1(i+1,j,n)     ; interpolate in x
               bydip=fx1*by1(i,j,n)+fx0*by1(i+1,j,n)
               bm0(kz)=sqrt(bxdip^2+bydip^2)
               dx(kz)=delx*bxdip/bm0(kz)
               dy(kz)=dely*bydip/bm0(kz)
               x1=x0(kz)-dx(kz)
               y1=y0(kz)-dy(kz)
               xx=(x1-xmin)/delx & i=long(xx)>0<(nx-1)
               fx0=xx-float(i) & fx1=1.0-fx0
               yy=(y1-ymin)/dely & j=long(yy)>0<(ny-1)
               fy0=yy-float(j) & fy1=1.0-fy0
               b1=fy1*(fx1*bz1(i  ,j  ,n)+fx0*bz1(i+1,j  ,n))  $
                 +fy0*(fx1*bz1(i  ,j+1,n)+fx0*bz1(i+1,j+1,n))
               x2=x0(kz)+dx(kz)
               y2=y0(kz)+dy(kz)
               xx=(x2-xmin)/delx & i=long(xx)>0<(nx-1)
               fx0=xx-float(i) & fx1=1.0-fx0
               yy=(y2-ymin)/dely & j=long(yy)>0<(ny-1)
               fy0=yy-float(j) & fy1=1.0-fy0
               b2=fy1*(fx1*bz1(i  ,j  ,n)+fx0*bz1(i+1,j  ,n))  $
                 +fy0*(fx1*bz1(i  ,j+1,n)+fx0*bz1(i+1,j+1,n))
               dbz(kz)=b2-b1
             endfor
           endif
           end
        1: begin  ; y-ribs
           index=where(bz1(i1:i2,j1:j2-1,n)*bz1(i1:i2,j1+1:j2,n) lt 0.0,nzero)
           if nzero gt 0 then begin
             x0=fltarr(nzero)  &  y0=fltarr(nzero)  &  z0=fltarr(nzero)
             dx=fltarr(nzero)  &  dy=fltarr(nzero)
             bm0=fltarr(nzero)
             dbz=fltarr(nzero)
             for kz=0,nzero-1 do begin
               i=i1+(index(kz) mod (i2-i1+1))
               j=j1+index(kz)/(i2-i1+1)
               fy0=bz1(i,j,n)/(bz1(i,j,n)-bz1(i,j+1,n))
               fy1=1.0-fy0
               x0(kz)=xmin+delz*float(i)
               y0(kz)=ymin+dely*(float(j)+fy0)
               z0(kz)=zmin+delz*float(n)
               bxdip=fy1*bx1(i,j,n)+fy0*bx1(i,j+1,n)     ; interpolate in y
               bydip=fy1*by1(i,j,n)+fy0*by1(i,j+1,n)
               bm0(kz)=sqrt(bxdip^2+bydip^2)
               dx(kz)=delx*bxdip/bm0(kz)
               dy(kz)=dely*bydip/bm0(kz)
               x1=x0(kz)-dx(kz)
               y1=y0(kz)-dy(kz)
               xx=(x1-xmin)/delx & i=long(xx)>0<(nx-1)
               fx0=xx-float(i) & fx1=1.0-fx0
               yy=(y1-ymin)/dely & j=long(yy)>0<(ny-1)
               fy0=yy-float(j) & fy1=1.0-fy0
               b1=fy1*(fx1*bz1(i  ,j  ,n)+fx0*bz1(i+1,j  ,n))  $
                 +fy0*(fx1*bz1(i  ,j+1,n)+fx0*bz1(i+1,j+1,n))
               x2=x0(kz)+dx(kz)
               y2=y0(kz)+dy(kz)
               xx=(x2-xmin)/delx & i=long(xx)>0<(nx-1)
               fx0=xx-float(i) & fx1=1.0-fx0
               yy=(y2-ymin)/dely & j=long(yy)>0<(ny-1)
               fy0=yy-float(j) & fy1=1.0-fy0
               b2=fy1*(fx1*bz1(i  ,j  ,n)+fx0*bz1(i+1,j  ,n))  $
                 +fy0*(fx1*bz1(i  ,j+1,n)+fx0*bz1(i+1,j+1,n))
               dbz(kz)=b2-b1
             endfor
           endif
           end
        2: begin  ; z-ribs
           index=where(bz1(i1:i2,j1:j2,n)*bz1(i1:i2,j1:j2,n+1) lt 0.0,nzero)
           if nzero gt 0 then begin
             x0=fltarr(nzero)  &  y0=fltarr(nzero)  &  z0=fltarr(nzero)
             dx=fltarr(nzero)  &  dy=fltarr(nzero)
             bm0=fltarr(nzero)
             dbz=fltarr(nzero)
             for kz=0,nzero-1 do begin
               i=i1+(index(kz) mod (i2-i1+1))
               j=j1+index(kz)/(i2-i1+1)
               fz0=bz1(i,j,n)/(bz1(i,j,n)-bz1(i,j,n+1))
               fz1=1.0-fz0
               x0(kz)=xmin+delx*float(i)
               y0(kz)=ymin+dely*float(j)
               z0(kz)=zmin+delz*(float(n)+fz0)
               bxdip=fz1*bx1(i,j,n)+fz0*bx1(i,j,n+1)     ; interpolate in z
               bydip=fz1*by1(i,j,n)+fz0*by1(i,j,n+1)
               bm0(kz)=sqrt(bxdip^2+bydip^2)
               dx(kz)=delx*bxdip/bm0(kz)
               dy(kz)=dely*bydip/bm0(kz)
               x1=x0(kz)-dx(kz)
               y1=y0(kz)-dy(kz)
               xx=(x1-xmin)/delx & i=long(xx)>0<(nx-1)
               fx0=xx-float(i) & fx1=1.0-fx0
               yy=(y1-ymin)/dely & j=long(yy)>0<(ny-1)
               fy0=yy-float(j) & fy1=1.0-fy0
               b1=fz1*(fy1*(fx1*bz1(i  ,j  ,n  )+fx0*bz1(i+1,j  ,n  ))  $
                      +fy0*(fx1*bz1(i  ,j+1,n  )+fx0*bz1(i+1,j+1,n  ))) $
                 +fz0*(fy1*(fx1*bz1(i  ,j  ,n+1)+fx0*bz1(i+1,j  ,n+1))  $
                      +fy0*(fx1*bz1(i  ,j+1,n+1)+fx0*bz1(i+1,j+1,n+1)))
               x2=x0(kz)+dx(kz)
               y2=y0(kz)+dy(kz)
               xx=(x2-xmin)/delx & i=long(xx)>0<(nx-1)
               fx0=xx-float(i) & fx1=1.0-fx0
               yy=(y2-ymin)/dely & j=long(yy)>0<(ny-1)
               fy0=yy-float(j) & fy1=1.0-fy0
               b2=fz1*(fy1*(fx1*bz1(i  ,j  ,n  )+fx0*bz1(i+1,j  ,n  ))  $
                      +fy0*(fx1*bz1(i  ,j+1,n  )+fx0*bz1(i+1,j+1,n  ))) $
                 +fz0*(fy1*(fx1*bz1(i  ,j  ,n+1)+fx0*bz1(i+1,j  ,n+1))  $
                      +fy0*(fx1*bz1(i  ,j+1,n+1)+fx0*bz1(i+1,j+1,n+1)))
               dbz(kz)=b2-b1
             endfor
           endif
           end
        endcase
        if nzero gt 0 then begin
          k=where(dbz gt amin*bm0,ndip)
          if ndip gt 0 then begin
            dx=0.5*dx(k)        ; reduce length for plotting bars
            dy=0.5*dy(k)
            if num1 eq 0 then begin
              x2a=x0(k)-dx  &  x2b=x0(k)+dx
              y2a=y0(k)-dy  &  y2b=y0(k)+dy
              z2a=z0(k)     &  z2b=z0(k)
              ccc=lonarr(ndip)+cc3
              ind=lindgen(ndip)
            endif else begin
              x2a=[x2a,x0(k)-dx]  &  x2b=[x2b,x0(k)+dx]
              y2a=[y2a,y0(k)-dy]  &  y2b=[y2b,y0(k)+dy]
              z2a=[z2a,z0(k)]     &  z2b=[z2b,z0(k)]
              ccc=[ccc,lonarr(ndip)+cc3]
              ind=[ind,lindgen(ndip)+num2]
            endelse
            num2=num2+ndip
            num1=num1+ndip
            ndips(m)=ndips(m)+ndip
          endif
        endif
      endfor    ; end loop over three cases (x,y,z)
      print,string(n,fblue,ndips,  $
            format='("  n=",i2,", fblue=",f4.2,", ndips=",3(i4,2x))')
      ndips_tot=ndips_tot+long(total(ndips))
      if ndips_tot eq 0 then nmx=(nmx-1)>(nzmin+2)
    endfor       ; end loop over heights
  endif
endif
print,string(num1,format='("Number of segments: num1=",i6)')
;
;  Eliminate "snipped sections", convert to viewing coordinates, and plot
;  sections in order of decreasing distance to observer:
;
if num1 gt 0 then begin
  subr_rotate,alt,azim,x2a(ind),y2a(ind),z2a(ind),xsa,ysa,zsa
  subr_rotate,alt,azim,x2b(ind),y2b(ind),z2b(ind),xsb,ysb,zsb
  ccc=ccc(ind)
  if num1 gt 0 then begin
    pos=0.5*(zsa+zsb)
    ind=sort(pos)
    if keyword_set(flbw) then begin
      for nn=0L,num1-1 do begin
        n=ind(nn)
        plots,[xsa(n),xsb(n)],[ysa(n),ysb(n)],/data,  $
              color=cc2,thick=th2bw,noclip=0
        plots,[xsa(n),xsb(n)],[ysa(n),ysb(n)],/data,  $
              color=cc1,thick=th1bw,noclip=0
      endfor
    endif else begin
      for nn=0L,num1-1 do begin
        n=ind(nn)
        plots,[xsa(n),xsb(n)],[ysa(n),ysb(n)],/data,  $
              color=ccc(n),thick=thcol,noclip=0
      endfor
    endelse
  endif
endif
;
;  Close Postcript file:
;
if keyword_set(epsfile) then begin
  device,/close
  print,'Plot saved on '+epsfile
  set_plot,'X'
  hexa_colour
endif
;
print,'-----'
end

;------------------------------------------------------------------
;  Convert "box" coordinates (x,y,z) into coordinates (xv,yv,zv)
;  using viewing angles (alt,azim) in degrees:
;------------------------------------------------------------------
pro subr_rotate,alt,azim,x,y,z,xv,yv,zv
@hexa.blk
;
azimr=azim/!radeg
xv =x*cos(azimr)-y*sin(azimr)
yv1=x*sin(azimr)+y*cos(azimr)
;
theta=(90.0-alt)/!radeg
zv=z*cos(theta)-yv1*sin(theta)
yv=z*sin(theta)+yv1*cos(theta)
;
end

;------------------------------------------------------------------
;  Convert viewed coordinates (xv,yv) into "box" coordinates (x,y) at z=0,
;  using viewing angles (alt,azim) in degrees:
;------------------------------------------------------------------
pro subr_invert1,alt,azim,xv,yv,x,y,ind,zlevel=zlevel
@hexa.blk
;
if n_elements(zlevel) eq 0 then zlevel=0.0
;
theta=(90.0-alt)/!radeg
yv1=(yv-zlevel*sin(theta))/cos(theta)
;
azimr=azim/!radeg
x= xv*cos(azimr)+yv1*sin(azimr)
y=-xv*sin(azimr)+yv1*cos(azimr)
;
ind=where(x gt xmin and x lt xmax and y gt ymin and y lt ymax)
;
end

;------------------------------------------------------------------
;  Convert heliocentric coordinates (xs,ys,zs) into "box" coordinates
;  (x,y,z) using viewing angles (lon,lat) in degrees:
;------------------------------------------------------------------
pro subr_invert2,alt,azim,xv,yv,zv,x,y,z,flg
;@hexa.blk
;
latr=lat/!radeg
xs1= xs
ys1= ys*cos(latr)+zs*sin(latr)
zs1=-ys*sin(latr)+zs*cos(latr)
;
lonr=lon/!radeg
tmp=-xs1*sin(lonr)+zs1*cos(lonr)
xs1= xs1*cos(lonr)+zs1*sin(lonr)
zs1=tmp
;
r=sqrt(xs1^2+ys1^2+zs1^2)
t=acos(ys1/r)
p=atan(xs1,zs1)
;
x=p/delt0
y=-alog(tan(0.5*t))/delt0
z=(alog(r)/delt0)>zmin<zmax
;
flg=(x ge xmin and x le xmax and y ge ymin and y le ymax)
end

;------------------------------------------------------------------
;  Show surface evolution:
;------------------------------------------------------------------
pro prog_evolve
@hexa.blk
;
answer=''
;read,prompt='Close XINTERANIMATE widget? (y/n) ',answer
;if answer eq 'y' then xinteranimate,/close
;
;  Open file with time-dependent boundary conditions:
;
file=dir+filename+'_evolve'
openr,1,file,/f77_unform
opt=0L
readu,1,opt
if opt ne 1L then begin
  print,'Error reading '+file
  close,1
  return
endif
print,'Open '+file
ax=fltarr(nx,ny+1)
ay=fltarr(nx+1,ny)
num=nmajor*nminor
bmax=0.0
for n=0,num-1 do begin
  readu,1,ax
  readu,1,ay
  bz=(ay(1:nx,*)-ay(0:nx-1,*))/delx-(ax(*,1:ny)-ax(*,0:ny-1))/dely
  bmax=max([bmax,max(abs(bz))])/2.0
endfor
bmin=-bmax
close,1
openr,1,file,/f77_unform
readu,1,opt
;
;  Create movie display:
;
zoom=2
window,4,xs=zoom*nx,ys=zoom*ny
xinteranimate,set=[zoom*nx,zoom*ny,num],/showload,/track, $
              title='Surface Evolution'

aax0=fltarr(nx,ny+1,num)
aay0=fltarr(nx+1,ny,num)

for n=0,num-1 do begin
  readu,1,ax
  aax0(*,*,n)=ax
  readu,1,ay
  aay0(*,*,n)=ay
  bz=(ay(1:nx,*)-ay(0:nx-1,*))/delx-(ax(*,1:ny)-ax(*,0:ny-1))/dely
  img=bytscl(bz,min=bmin,max=bmax)
  if zoom gt 1 then img=rebin(img,zoom*nx,zoom*ny,/sample)
  wset,4
  tv,img
  xinteranimate,frame=n,window=4
endfor
close,1
;
xinteranimate,30.
end

;------------------------------------------------------------------
;  Reset field-line counter:
;------------------------------------------------------------------
pro subr_reset
@hexa.trc
;
np=0
xp=fltarr(1)
yp=fltarr(1)
zp=fltarr(1)
;
;nlos=0
;ip1=intarr(1)
;ip2=intarr(1)
;
end

;------------------------------------------------------------------
;  Trace field line from starting point (x0,y0,z0):
;------------------------------------------------------------------
pro trace,x0,y0,z0,zbase=zbase,smin=smin,smax=smax
@hexa.blk
@hexa.trc
;
;  Minimum and maximum position along field line:
;
if n_elements(smin) eq 0 then smin=-999.
smin=smin<0.
if n_elements(smax) eq 0 then smax= 999.
smax=smax>0.
if smin eq smax then begin
  print,'Warning: smin=smax=0, set to [-0.5,+0.5]'
  smin=-0.5  &  smax=0.5
endif
;
;  Lower boundary:
;
if n_elements(zbase) eq 0 then zbase=zmin
;
;  Check that starting point lies within the box or on lower boundary:
;
dm=0.1*delx   ; minimum distance from front, back and top boundaries
if x0 gt xmax    then x0=xmax
if x0 lt xmin    then x0=xmin
if y0 gt ymax-dm then y0=ymax-dm
if y0 lt ymin+dm then y0=ymin+dm
if z0 gt zmax-dm then z0=zmax-dm
if z0 lt zbase   then z0=zbase
;
;  Initialize arrays:
;
kdim=100
xtrc=fltarr(kdim)
ytrc=fltarr(kdim)
ztrc=fltarr(kdim)
strc=fltarr(kdim)
ds2=fltarr(kdim)
bx2=fltarr(kdim)
by2=fltarr(kdim)
bz2=fltarr(kdim)
bm2=fltarr(kdim)
dum=fltarr(100)
;
;  Trace a field line starting from an arbitrary point, (x0,y0,z0).
;  On exit, kp is the index of the starting point.
;
kp=kdim-1
k=kp
xtrc(k)=x0
ytrc(k)=y0
ztrc(k)=z0
strc(k)=0.0
tr_magn,xtrc(k),ytrc(k),ztrc(k),ds2a,bx2a,by2a,bz2a
ds2(k)=ds2a
bx2(k)=bx2a
by2(k)=by2a
bz2(k)=bz2a
bm2(k)=sqrt(bx2(k)^2+by2(k)^2+bz2(k)^2)
;
;  Trace backward:
;
if (ztrc(k) gt zbase) or (ztrc(k) eq zbase and bz2(k) lt 0.0) then begin
  ds2(k-1)=ds2(k)
  mode=0
  n=0
  repeat begin
    n=n+1
    k=k-1
;
;  If necessary, increase array size (backward):
;
    if k eq 0 then begin
      xtrc=[dum,xtrc]
      ytrc=[dum,ytrc]
      ztrc=[dum,ztrc]
      strc=[dum,strc]
      ds2=[dum,ds2]
      bx2=[dum,bx2]
      by2=[dum,by2]
      bz2=[dum,bz2]
      bm2=[dum,bm2]
      kdim=kdim+100
      kp=kp+100
      k=k+100
    endif
    k1=k+1
;
;  Compute magnetic field at half-step (backward):
;
    fact=0.5*ds2(k)/bm2(k1)
    x2h=xtrc(k1)-fact*bx2(k1)
    y2h=ytrc(k1)-fact*by2(k1)
    z2h=ztrc(k1)-fact*bz2(k1)
    tr_magn,x2h,y2h,z2h,ds2h,bx2h,by2h,bz2h
    bm2h=sqrt(bx2h^2+by2h^2+bz2h^2)
;
;  On the first backward step, if the field line is curved downward,
;  make sure that the step size is smaller than about one fifth
;  of the (estimated) loop length, otherwise reduce step size and
;  recompute the magnetic field at the half-step:
;
    if n eq 1 then begin
      e1=bz2(k1)/bm2(k1)
      eh=bz2h   /bm2h
      deds=(e1-eh)/(0.5*ds2(k))
      if deds lt 0.0 then begin
;        rad=exp(delt0*ztrc(k1))
;        length=-2.0/deds*sqrt(e1^2-2.*rad*deds)
        length=-2.0/deds*sqrt(e1^2-2.*deds)
        if ds2(k) gt 0.21*length then begin
          ds2(k)=0.21*length
          mode=1
          fact=0.5*ds2(k)/bm2(k1)
          x2h=xtrc(k1)-fact*bx2(k1)
          y2h=ytrc(k1)-fact*by2(k1)
          z2h=ztrc(k1)-fact*bz2(k1)
          tr_magn,x2h,y2h,z2h,ds2h,bx2h,by2h,bz2h
          bm2h=sqrt(bx2h^2+by2h^2+bz2h^2)
        endif
      endif
    endif
;
;  Compute magnetic field at full-step (backward):
;
    fact=ds2(k)/bm2h
    xtrc(k)=xtrc(k1)-fact*bx2h
    ytrc(k)=ytrc(k1)-fact*by2h
    ztrc(k)=ztrc(k1)-fact*bz2h
    tr_magn,xtrc(k),ytrc(k),ztrc(k),ds2a,bx2a,by2a,bz2a
    ds2(k-1)=ds2a
    bx2(k)=bx2a
    by2(k)=by2a
    bz2(k)=bz2a
    bm2(k)=sqrt(bx2(k)^2+by2(k)^2+bz2(k)^2)
    strc(k)=strc(k1)-ds2(k)
    if mode then ds2(k-1)=ds2(k)
  endrep until ((strc(k) lt smin)  or (strc(k) gt smax) or  $
                (ztrc(k) le zbase) or (ztrc(k) ge zmax) or  $
                (ytrc(k) le ymin)  or (ytrc(k) ge ymax) or  $
                (xtrc(k) lt x2min) or (xtrc(k) gt x2max) or (N_elements(xtrc) gt 10000l))
;
;  Correct last point of backward trace:
;
  if ztrc(k) lt zbase then begin
    frac=(zbase-ztrc(k1))/(ztrc(k)-ztrc(k1))
    frac1=1.0-frac
    xtrc(k)=xtrc(k1)*frac1+xtrc(k)*frac
    ytrc(k)=ytrc(k1)*frac1+ytrc(k)*frac
    ztrc(k)=ztrc(k1)*frac1+ztrc(k)*frac
    strc(k)=strc(k1)*frac1+strc(k)*frac
    ds2(k)=frac*ds2(k)
    bx2(k)=bx2(k1)*frac1+bx2(k)*frac
    by2(k)=by2(k1)*frac1+by2(k)*frac
    bz2(k)=bz2(k1)*frac1+bz2(k)*frac
    bm2(k)=bm2(k1)*frac1+bm2(k)*frac
  endif
endif
;
;  Reorganize data arrays:
;
k0=k
xtrc=[xtrc(k0:kdim-1),dum]
ytrc=[ytrc(k0:kdim-1),dum]
ztrc=[ztrc(k0:kdim-1),dum]
strc=[strc(k0:kdim-1),dum]
ds2=[ds2(k0:kdim-1),dum]
bx2=[bx2(k0:kdim-1),dum]
by2=[by2(k0:kdim-1),dum]
bz2=[bz2(k0:kdim-1),dum]
bm2=[bm2(k0:kdim-1),dum]
kdim=kdim-k0+100
kp=kp-k0
k=kp
;
;  Trace forward:
;
if (ztrc(k) gt zbase) or (ztrc(k) eq zbase and bz2(k) gt 0.0) then begin
  mode=0
  n=0
  repeat begin
    n=n+1
    k=k+1
;
;  If necessary, increase array size (forward):
;
    if k eq kdim-1 then begin
      xtrc=[xtrc,dum]
      ytrc=[ytrc,dum]
      ztrc=[ztrc,dum]
      strc=[strc,dum]
      ds2=[ds2,dum]
      bx2=[bx2,dum]
      by2=[by2,dum]
      bz2=[bz2,dum]
      bm2=[bm2,dum]
      kdim=kdim+100
    endif
    k1=k-1
;
;  Compute magnetic field at half-step (forward):
;
    fact=0.5*ds2(k1)/bm2(k1)
    x2h=xtrc(k1)+fact*bx2(k1)
    y2h=ytrc(k1)+fact*by2(k1)
    z2h=ztrc(k1)+fact*bz2(k1)
    tr_magn,x2h,y2h,z2h,ds2h,bx2h,by2h,bz2h
    bm2h=sqrt(bx2h^2+by2h^2+bz2h^2)
;
;  On the first forward step, if the field line is curved downward,
;  make sure that the step size is smaller than about one fifth
;  of the (estimated) loop length, otherwise reduce step size and
;  recompute the magnetic field at the half-step:
;
    if n eq 1 then begin
      e1=bz2(k1)/bm2(k1)
      eh=bz2h   /bm2h
      deds=(eh-e1)/(0.5*ds2(k1))
      if deds lt 0.0 then begin
;        rad=exp(delt0*ztrc(k1))
;        length=-2.0/deds*sqrt(e1^2-2.*rad*deds)
        length=-2.0/deds*sqrt(e1^2-2.*deds)
        if ds2(k1) gt 0.21*length then begin
          ds2(k1)=0.21*length
          mode=1
          fact=0.5*ds2(k1)/bm2(k1)
          x2h=xtrc(k1)+fact*bx2(k1)
          y2h=ytrc(k1)+fact*by2(k1)
          z2h=ztrc(k1)+fact*bz2(k1)
          tr_magn,x2h,y2h,z2h,ds2h,bx2h,by2h,bz2h
          bm2h=sqrt(bx2h^2+by2h^2+bz2h^2)
        endif
      endif
    endif
;
;  Compute magnetic field at full-step (backward):
;
    fact=ds2(k1)/bm2h
    xtrc(k)=xtrc(k1)+fact*bx2h
    ytrc(k)=ytrc(k1)+fact*by2h
    ztrc(k)=ztrc(k1)+fact*bz2h
    tr_magn,xtrc(k),ytrc(k),ztrc(k),ds2a,bx2a,by2a,bz2a
    ds2(k)=ds2a
    bx2(k)=bx2a
    by2(k)=by2a
    bz2(k)=bz2a
    bm2(k)=sqrt(bx2(k)^2+by2(k)^2+bz2(k)^2)
    strc(k)=strc(k1)+ds2(k1)
    if mode then ds2(k)=ds2(k1)
  endrep until ((strc(k) lt smin)  or (strc(k) gt smax) or  $
                (ztrc(k) le zbase) or (ztrc(k) ge zmax) or  $
                (ytrc(k) le ymin)  or (ytrc(k) ge ymax) or  $
                (xtrc(k) lt x2min) or (xtrc(k) gt x2max) or (N_elements(xtrc) gt 10000l))
;
;  Correct last point of backward trace:
;
  if ztrc(k) le zbase then begin
    frac=(zbase-ztrc(k1))/(ztrc(k)-ztrc(k1))
    frac1=1.0-frac
    xtrc(k)=xtrc(k1)*frac1+xtrc(k)*frac
    ytrc(k)=ytrc(k1)*frac1+ytrc(k)*frac
    ztrc(k)=ztrc(k1)*frac1+ztrc(k)*frac
    strc(k)=strc(k1)*frac1+strc(k)*frac
    ds2(k)=frac*ds2(k)
    bx2(k)=bx2(k1)*frac1+bx2(k)*frac
    by2(k)=by2(k1)*frac1+by2(k)*frac
    bz2(k)=bz2(k1)*frac1+bz2(k)*frac
    bm2(k)=bm2(k1)*frac1+bm2(k)*frac
  endif
endif
;
;  Reorganize data arrays:
;
xtrc=xtrc(0:k)
ytrc=ytrc(0:k)
ztrc=ztrc(0:k)
strc=strc(0:k)
ds2=ds2(0:k)
bx2=bx2(0:k)
by2=by2(0:k)
bz2=bz2(0:k)
bm2=bm2(0:k)
;
end

;----------------------------------------------------------------
;  Linear interpolation in 3D model:
;----------------------------------------------------------------
pro tr_magn,xx2,yy2,zz2,ds2,fx2,fy2,fz2,par=par,verbose=verbose
@hexa.blk
;
x2=xx2
xdel=xmax-xmin
while (x2 lt xmin) do x2=x2+xdel
while (x2 gt xmax) do x2=x2-xdel
dx2=(x2-xmin)/delx
i=fix(dx2)<(nx-1)  &  ip=i+1
fx0=dx2-float(i)
fx1=1.0-fx0
;
y2=yy2
dy2=(y2-ymin)/dely
j=fix(dy2)<(ny-1)  &  jp=j+1
fy0=dy2-float(j)
fy1=1.0-fy0
;
z2=zz2>0.
dz2=(z2-zmin)/delz
k=fix(dz2)>0<(nz-1)  &  kp=k+1
fz0=(z2-k*delz)/delz
fz1=1.0-fz0
;
ds2=0.3*delx
;
;  Parameter:  0=MAG, 1=CUR, 2=VEL, 3=FOR, 4=alpha, 5=grad, 6=div4
;
if n_elements(par) eq 0 then par=0
case par of
0: begin
   fx2=fz1*(fy1*(fx1*bx1(i ,j ,k )+fx0*bx1(ip,j ,k ))    $
           +fy0*(fx1*bx1(i ,jp,k )+fx0*bx1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*bx1(i ,j ,kp)+fx0*bx1(ip,j ,kp))    $
           +fy0*(fx1*bx1(i ,jp,kp)+fx0*bx1(ip,jp,kp)))
   fy2=fz1*(fy1*(fx1*by1(i ,j ,k )+fx0*by1(ip,j ,k ))    $
           +fy0*(fx1*by1(i ,jp,k )+fx0*by1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*by1(i ,j ,kp)+fx0*by1(ip,j ,kp))    $
           +fy0*(fx1*by1(i ,jp,kp)+fx0*by1(ip,jp,kp)))
   fz2=fz1*(fy1*(fx1*bz1(i ,j ,k )+fx0*bz1(ip,j ,k ))    $
           +fy0*(fx1*bz1(i ,jp,k )+fx0*bz1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*bz1(i ,j ,kp)+fx0*bz1(ip,j ,kp))    $
           +fy0*(fx1*bz1(i ,jp,kp)+fx0*bz1(ip,jp,kp)))
   end
1: begin
   fx2=fz1*(fy1*(fx1*cx1(i ,j ,k )+fx0*cx1(ip,j ,k ))    $
           +fy0*(fx1*cx1(i ,jp,k )+fx0*cx1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*cx1(i ,j ,kp)+fx0*cx1(ip,j ,kp))    $
           +fy0*(fx1*cx1(i ,jp,kp)+fx0*cx1(ip,jp,kp)))
   fy2=fz1*(fy1*(fx1*cy1(i ,j ,k )+fx0*cy1(ip,j ,k ))    $
           +fy0*(fx1*cy1(i ,jp,k )+fx0*cy1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*cy1(i ,j ,kp)+fx0*cy1(ip,j ,kp))    $
           +fy0*(fx1*cy1(i ,jp,kp)+fx0*cy1(ip,jp,kp)))
   fz2=fz1*(fy1*(fx1*cz1(i ,j ,k )+fx0*cz1(ip,j ,k ))    $
           +fy0*(fx1*cz1(i ,jp,k )+fx0*cz1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*cz1(i ,j ,kp)+fx0*cz1(ip,j ,kp))    $
           +fy0*(fx1*cz1(i ,jp,kp)+fx0*cz1(ip,jp,kp)))
   end
2: begin
   fx2=fz1*(fy1*(fx1*vx1(i ,j ,k )+fx0*vx1(ip,j ,k ))    $
           +fy0*(fx1*vx1(i ,jp,k )+fx0*vx1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*vx1(i ,j ,kp)+fx0*vx1(ip,j ,kp))    $
           +fy0*(fx1*vx1(i ,jp,kp)+fx0*vx1(ip,jp,kp)))
   fy2=fz1*(fy1*(fx1*vy1(i ,j ,k )+fx0*vy1(ip,j ,k ))    $
           +fy0*(fx1*vy1(i ,jp,k )+fx0*vy1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*vy1(i ,j ,kp)+fx0*vy1(ip,j ,kp))    $
           +fy0*(fx1*vy1(i ,jp,kp)+fx0*vy1(ip,jp,kp)))
   fz2=fz1*(fy1*(fx1*vz1(i ,j ,k )+fx0*vz1(ip,j ,k ))    $
           +fy0*(fx1*vz1(i ,jp,k )+fx0*vz1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*vz1(i ,j ,kp)+fx0*vz1(ip,j ,kp))    $
           +fy0*(fx1*vz1(i ,jp,kp)+fx0*vz1(ip,jp,kp)))
   end
3: begin
   fx2=fz1*(fy1*(fx1*ffx1(i ,j ,k )+fx0*ffx1(ip,j ,k ))    $
           +fy0*(fx1*ffx1(i ,jp,k )+fx0*ffx1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*ffx1(i ,j ,kp)+fx0*ffx1(ip,j ,kp))    $
           +fy0*(fx1*ffx1(i ,jp,kp)+fx0*ffx1(ip,jp,kp)))
   fy2=fz1*(fy1*(fx1*ffy1(i ,j ,k )+fx0*ffy1(ip,j ,k ))    $
           +fy0*(fx1*ffy1(i ,jp,k )+fx0*ffy1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*ffy1(i ,j ,kp)+fx0*ffy1(ip,j ,kp))    $
           +fy0*(fx1*ffy1(i ,jp,kp)+fx0*ffy1(ip,jp,kp)))
   fz2=fz1*(fy1*(fx1*ffz1(i ,j ,k )+fx0*ffz1(ip,j ,k ))    $
           +fy0*(fx1*ffz1(i ,jp,k )+fx0*ffz1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*ffz1(i ,j ,kp)+fx0*ffz1(ip,j ,kp))    $
           +fy0*(fx1*ffz1(i ,jp,kp)+fx0*ffz1(ip,jp,kp)))
   end
4: begin
   fx2=fz1*(fy1*(fx1*alpha1(i ,j ,k )+fx0*alpha1(ip,j ,k ))    $
           +fy0*(fx1*alpha1(i ,jp,k )+fx0*alpha1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*alpha1(i ,j ,kp)+fx0*alpha1(ip,j ,kp))    $
           +fy0*(fx1*alpha1(i ,jp,kp)+fx0*alpha1(ip,jp,kp)))
   fy2=0.
   fz2=0.
   end
5: begin
   fx2=fz1*(fy1*(fx1*grad1(i ,j ,k )+fx0*grad1(ip,j ,k ))    $
           +fy0*(fx1*grad1(i ,jp,k )+fx0*grad1(ip,jp,k )))   $
      +fz0*(fy1*(fx1*grad1(i ,j ,kp)+fx0*grad1(ip,j ,kp))    $
           +fy0*(fx1*grad1(i ,jp,kp)+fx0*grad1(ip,jp,kp)))
   fy2=0.
   fz2=0.
   end
6: begin
   fx2=fz1*(fy1*(fx1*div4(i ,j ,k )+fx0*div4(ip,j ,k ))    $
           +fy0*(fx1*div4(i ,jp,k )+fx0*div4(ip,jp,k )))   $
      +fz0*(fy1*(fx1*div4(i ,j ,kp)+fx0*div4(ip,j ,kp))    $
           +fy0*(fx1*div4(i ,jp,kp)+fx0*div4(ip,jp,kp)))
   fy2=0.
   fz2=0.
   end
endcase
;
if keyword_set(verbose) then begin
  var=['MAG','CUR','VEL','FOR','ALPHA','GRAD','DIV4']
  print,string(x2,y2,z2,  $
        format='("x2=",f7.2,", y2=",f7.2,", z2=",f7.2)')
  print,string(k,nx,ny,  $
        format='("k=",i3,", nx=",i4,", ny=",i4)')
  print,string(i,j,format='("i=",i4,", j=",i4)')
  print,string(fx0,fy0,fz0,  $
        format='("fx0=",f7.3,", fy0=",f7.3,", fz0=",f7.3)')
  print,'var=',var(par)
  print,string(fx2,fy2,fz2,  $
        format='("fx2=",e10.3,", fy2=",e10.3,", fz2=",e10.3)')
endif
end

;------------------------------------------------------------------
;  Save field-line starting points:
;------------------------------------------------------------------
pro prog_save_lines
@hexa.blk
@hexa.trc
;
if n_elements(np) eq 0 then np=0
if np eq 0 then begin
  print,'No field lines to save'
  print,'-----'
  return
endif
answer=''
read,prompt='Enter field line ID (0=quit):  ',answer
if answer eq '' or answer eq '0' then begin
  print,'Quit'
  print,'-----'
  return
endif
flid=fix(answer)>1
file=dir+string(filename,nt,format='(a,"_",i5.5,"_P")')  $
        +strcompress(flid,/remove_all)
dum=findfile(file,count=count)
if count gt 0 then begin
  print,'Warning: file '+file+' already exists'
  answer=''
  read,prompt='Overwrite file? (y/n,def=n)  ',answer
  if answer ne 'y' then begin
    print,'Not saved'
    print,'-----'
    return
  endif
endif
print,'Save data in '+file
openw,1,file
printf,1,np
for n=0,np-1 do printf,1,string(xp(n),yp(n),zp(n),format='(3(f8.3,1x))')
close,1
print,'-----'
end














function add_element,x,val,n
if n eq 1 then begin
xx = intarr(n)
xx(n-1) = val
endif else begin
xx=intarr(n)
xx(0:n-2) = x
xx(n-1) = val
endelse

return,xx
end




function tension,bbz1
@hexa.blk
tz=fltarr(nx+1,ny+1,nz+1)
dxb = (bbz1(1:nx,*,*)-bbz1(0:nx-1,*,*))/delx
tz(1:nx-1,*,*) = tz(1:nx-1,*,*) + bx1(1:nx-1,*,*)*( dxb(1:nx-1,*,*)+dxb(0:nx-2,*,*) )/2.;

dyb = (bbz1(*,1:ny,*)-bbz1(*,0:ny-1,*))/dely
tz(*,1:ny-1,*) = tz(*,1:ny-1,*) + by1(*,1:ny-1,*)*( dyb(*,1:ny-1,*) + dyb(*,0:ny-2,*))/2.

dzb = (bbz1(*,*,1:nz)-bbz1(*,*,0:nz-1))/delz
tz(*,*,1:nz-1) = tz(*,*,1:nz-1) + bz1(*,*,1:nz-1)*(dzb(*,*,1:nz-1)+dzb(*,*,0:nz-2))/2.


return,tz
end










pro frfinder,threshold=threshold;,3d=3d
@hexa.blk
@hexa.col
@hexa.trc
; magnetic pressure gradient

print,'Calculating magnetic pressure...'

pz = fltarr(nx+1,ny+1,nz+1)
px = pz
py = pz


b2 = bx1^2 + by1^2 + bz1^2
;Pz
pz2 = (b2(*,*,1:nz) - b2(*,*,0:nz-1))/delz ; grad (B^2) (z-ribs)
pz(*,*,1:nz-1) = (pz2(*,*,1:nz-1) + pz2(*,*,0:nz-2))/(-2.) ; grad(B^2) (avgd to corners)
pz = pz/2.


;Py
py2 =( b2(*,1:ny,*) - b2(*,0:ny-1,*))/dely
py(*,1:ny-1,0) = (py2(*,1:ny-1,*) + py2(*,0:ny-2,0))/(-2.)
py = py/2.

;Px
px2 =( b2(1:nx,*,*) - b2(0:nx-1,*,*))/delx
px(1:nx-1,*,*) = (px2(1:nx-1,*,*) + px2(0:nx-2,*,*))/(-2.)
px = px/2.



; magnetic tension
print,'Calculating magnetic tension...'

tz=tension(bz1)


tx=tension(bx1)

ty = tension(by1)



; seek flux rope:

print,'Looking for flux ropes...'
n=0

interval = 1
if not keyword_set(threshold) then threshold = 50.
bmax = threshold^2



xxx = intarr(1)
yyy = intarr(1)
zzz = intarr(1)
 for k=1,nz-2,interval do begin
  for j=1,ny-2,interval  do begin
   for i=1,nx-2,interval do begin

      if b2(i,j,k) ge bmax then begin ; if magnetic field is high...
         if tz(i,j,k) gt 0. then begin
           if tz(i,j,k+interval) lt 0 then begin
             if pz(i,j,k) lt 0 and pz(i,j,k+interval) gt 0 then begin
               if (px(i,j,k) lt 0 and px(i+interval,j,k) gt 0) or (py(i,j,k) lt 0 and py(i,j+interval,k) gt 0) then begin
                 if (tx(i,j,k) gt 0 and tx(i+interval,j,k) lt 0) or (ty(i,j,k) gt 0 and ty(i,j+interval,k) lt 0) then begin
                    n=n+1
                    xxx=add_element(xxx,i+0.5,n)
                    yyy=add_element(yyy,j+0.5,n)
                    zzz=add_element(zzz,k+0.5,n)
                 endif
               endif
             endif
           endif
         endif
       endif


;         interval=2
;         if b2(i,j,k) ge bmax then begin ; if magnetic field is high...
;          if tz(i,j,k-1) gt 0. and  tz(i,j,k+1) lt 0 then begin
;              if pz(i,j,k-1) lt 0 and pz(i,j,k+1) gt 0 then begin
;                if (px(i-1,j,k) lt 0 and px(i+1,j,k) gt 0) or (py(i,j-1,k) lt 0 and py(i,j+1,k) gt 0) then begin
;                  if (tx(i-1,j,k) gt 0 and tx(i+1,j,k) lt 0) or (ty(i,j-1,k) gt 0 and ty(i,j+1,k) lt 0) then begin
;                     n=n+1
;                     xxx=add_element(xxx,i,n)
;                     yyy=add_element(yyy,j,n)
;                     zzz=add_element(zzz,k,n)
;                  endif
;                endif
;              endif
;          endif
;        endif
    endfor
   endfor
  endfor
   print,'found ',n,' points belonging to flux rope axes!'
print,'maximum height of FR:',max(zzz),' pixels'
print,'maximum field strength probed:', threshold,'G'

nplus=0
nminus=0

xxp=intarr(1)
yyp=xxp
zzp=xxp
xxn=xxp
yyn=xxp
zzn=xxp

for l=0,n-1 do begin ; check alpha of FL
  i = xxx(l)
  j = yyy(l)
  k = zzz(l)

  alpha = cx1(i,j,k)*bx1(i,j,k) + cy1(i,j,k)*by1(i,j,k) + cz1(i,j,k)*bz1(i,j,k)
  alpha = alpha / (sqrt(b2(i,j,k)))
  alpha = alpha / sqrt(cx1(i,j,k)^2 + cy1(i,j,k)^2 + cz1(i,j,k)^2)

  if alpha gt 0 then begin
    nplus = nplus +1
    xxp = add_element(xxp,i,nplus)
    yyp = add_element(yyp,j,nplus)
    zzp = add_element(zzp,k,nplus)
  endif else begin
    nminus = nminus+1
    xxn = add_element(xxn,i,nminus)
    yyn = add_element(yyn,j,nminus)
    zzn = add_element(zzn,k,nminus)
  endelse
endfor

xp = float(xxx)*delx;*6./256.
yp = float(yyy)*dely;*6./256.
zp = float(zzz)*delz;*6./256.
np = n


window,5,xsize=500,ysize=500

;xx=findgen(nx+1)
;yy=findgen(ny+1)
;zz= findgen(nz+1)

contour,bz1(*,*,0),levels=[-300,-200,-100,-50,50,100,200,300],/isotropic,/nodata

  contour,bz1(*,*,0),levels=[-411,-205,-102,-51,-25.6],/overplot,/isotropic,color=green
  contour,bz1(*,*,0),levels=[25.6,51,102,205,411],/overplot,/isotropic,color=red

plots,xxp,yyp,psym=2,color=yellow ; positive alpha
plots,xxn,yyn,psym=2,color=magenta ; negative alpha


;if keyword_set(3d) then begin
window,6;,xsize=800,ysize=800
surface,dist(256),xrange=[0,nx],yrange=[0,ny],zrange=[0,nz],/nodata,xtitle='x',ytitle='y',ztitle='z',charsize=2,/save







plots,xxp,yyp,zzp,psym=2,/t3d,color=yellow,symsize=2
plots,xxn,yyn,zzn,psym=2,/t3d,color=magenta,symsize=2

for n=0,np-1 do begin
plots,[xxx(n),xxx(n)],[yyy(n),yyy(n)],[0,zzz(n)],/t3d,linestyle=0
endfor
;endif
contour,bz1(*,*,0),levels=[-300,-200,-100,-50,50,100,200,300],/nodata,/t3d,/noerase,zvalue=0,charsize=2

  contour,bz1(*,*,0),levels=[-411,-205,-102,-51,-25.6],color=green,/t3d,zvalue=0,/overplot
  contour,bz1(*,*,0),levels=[25.6,51,102,205,411],color=red,/t3d,zvalue=0,/overplot

end

;-----------------------------------------------------------------
;  Main level:
;-----------------------------------------------------------------
@hexa.blk
@hexa.col
@hexa.trc
;@hexa.qsl
;
common wid,wid_main
;
device,retain=2
hexa_colour        ; load color table
version='2.0'     ; version of HEXA software
cd,current=home   ; current directory

end
