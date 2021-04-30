pro preparation
@cleanup.blk

res=1

;########## Select size of simulation box ##############

print,'Magnetograms are ',strtrim(nx,2),'x',strtrim(ny,2),' pixels.'
jump4: print,'What resolution would you like the simulation to be run at?'
print,'(1) 128^3'
print,'(2) 256^3'
print,'(3) 512^3'
read,res

scale='n'

;sets resolution of box according to the choice made above
if res eq 1 then size=128 else if res eq 2 then size=256 else if res eq 3 then size=512 else begin
    print,'Invalid option! Please try again...'
    goto,jump4
endelse

;checks if box size is less than the size of the magnetogram
if max([nx,ny]) gt size then begin
        jump5: print,'WARNING: Magnetograms are larger than box size!'
        print,'Do you wish to continue? (y/n)'
        answer='a'
        read,answer
        if answer eq 'n' then begin
           goto,jump4
        endif else if answer ne 'y' then begin
           print, 'Invalid response! Please try again...'
           goto,jump5
        endif
endif
choice=0
print,''

;preview embedding magnetograms 
!p.background=255
window,0,xs=600,ys=300
f1=dist(301)
f2=fltarr(301,301)+1*max(f1)/2.
f2(50:249,50:249)=dist(200)*max(f1)/max(dist(200))
!p.multi=[0,2,1]
xx=findgen(301)/300.*6.
pg_plotimage,bytscl(f1),xx,xx,title='Resized',color=0,xrange=[0,6],yrange=[0,6],/isotropic
pg_plotimage,bytscl(f2),xx,xx,title='Embedded',color=0,xrange=[0,6],yrange=[0,6],/isotropic
!p.background=0


print,'--------------------------------------------------------------'
print,'Next, choose whether to resize the magnetograms to fit the'
print,'computational box, or to embed the magnetograms in the box'
print,'with a custom scaling. Examples of the two options are '
print,'displayed in Window 0.'
jump6: print, 'Do you want to'
print, '(1) Resize magnetogram to fill box'
print, '(2) Embed magnetogram in centre of box (custom scaling)'
read,choice

;Places magnetogram time series into box and rescales as appropriate
if choice eq 1 then begin
    resize,size,/resize ;resizes magnetogram to fit in box (see resize.pro)
endif else if choice eq 2 then begin
    resize,size,/embed ;allows the user to choose how big the magnetogram is to be in the box (see resize.pro)
endif else begin
    print,'Invalid option! Please try again...'
    goto,jump6
endelse


;displays a magnetogram inside the computational box so the user can decide if they are happy with their choices. If not, they can go back and change them.
answer='a'
print,'Window 0 previews frame '+strtrim(nt/2,2)+' within the box.'
jump8: print,'Do you wish redo the grid placement? (y/n)'
read,answer
if answer eq 'y' then begin
  goto,jump4
endif else if answer ne 'n' then begin
  print,'Invalid response! Please try again'
  goto,jump8
endif


print,'--------------------------------------------------------------'


answer='a'
print,'Do you want the computational box to be periodic? (y/n[recommended])'
read,answer
if answer eq 'y' then periodic=1 else periodic=0


answer='a'
print,'Do you want to have an open top boundary? (y/n[recommended])'
read,answer
if answer eq 'y' then open=1 else open=0


if open eq 1 then begin
  answer ='a'
  print,'Would you like to flux balance the magnetograms? (y/n)'
  read,answer
  cor_mag=mag_data
  choice = 0
  if answer eq 'n' then goto,fbskip
endif


;################### Choose and apply flux balancing method #####################
signtest=0
choice=0

;preview flux balancing
!p.background=255
!p.multi=0
window,0,xs=300,ys=300
centre=dist(192,192)
zro=max(centre)/2.
image=fltarr(256,256)+1*zro
image(8:247,8:247) = max(centre)
image(8+16:247-16,8+16:247-16) = zro
image(32:255-32,32:255-32) = centre
xx=findgen(256)/256.*6.

pg_plotimage,bytscl(image),xx,xx,/isotropic,color=0,xrange=[0,6],yrange=[0,6],title='Example of halo'
!p.background=0

if open eq 0 then begin
  print,'Now flux balancing must be applied to ensure that the net flux '
  print,'through the base is zero. '
endif


print,'Flux balancing can be achieved through either'
print,'adjusting the value of every non-zero valued pixel to achieve'
print,'flux balance, or through a halo of pixels around the edge of'
print,'the box (requires blank space around the magnetogram). An'
print,'example of the halo is displayed in Window 0.'
jump9: print,'Please choose the flux balancing method'
print,'(1) Adjust pixels'
print,'(2) Add halo around edge' 
read,choice


imb=fltarr(nt) ;flux imbalance for each frame
cor_mag=fltarr(size,size,nt) ;array to contain the flux-balanced magnetograms
perpix=fltarr(nt)

if choice eq 1 then begin ;adjust all pixels
   for t=0,nt-1 do begin
      dummy=mag_data(*,*,t)
      imb(t)=total(dummy);determine flux imbalance
      nonzero = where(abs(dummy) gt 25);find non-zero pixels
      n=n_elements(nonzero) ;count non-zero pixels
      perpix(t)=imb(t)/n ;flux imbalance per pixel
      if abs(perpix(t)) gt 25 then begin ;display error if some pixels are to change sign
        print,'FRAME ',strtrim(t,2),' WARNING: Some pixels change sign!'
        ;print,'Correction of ',strtrim(perpix,2),'G per pixel'
        signtest=signtest+1
      endif
      dummy(nonzero)=dummy(nonzero)-perpix(t) ;balance frame
      cor_mag(*,*,t)=dummy ; put frame into cor_mag
    endfor
    
    
    
    
    ;if the flux balancing changes the sign of some pixels, ask user if they want to continue or go back
    if signtest gt 0 then begin 
    window,0,xs=800,ys=500
    plot,perpix,xtitle='Frame',ytitle='Correction per pixel (G)',yrange=[min([-30,min(perpix)]),max([30,max(perpix)])],charsize=2
    oplot,[0,nt],[-25,-25],linestyle=2
    oplot,[0,nt],[25,25],linestyle=2
      jump10: Print,'Flux balancing caused some pixels to change sign in ',strtrim(signtest,2),' frames.'
      choice='a'
      print,'Do you wish to continue? (y/n)'
      read,choice
      wdelete,0
      if choice eq 'n' then goto,jump9 else if choice ne 'y' then begin
        print, 'Invalid response! Please try again...'
        goto,jump10
      endif
    endif
    
       
endif else if choice eq 2 then begin ;use halo to achieve flux balance

    width=size/16 ; width of halo (8 cells for 128, 16 for 256 etc...)
    ;The halo also has width/2 blank cells either side of it, so the
    ;full width of the halo is width*2 pixels

    if (xmargin lt width*4) or (ymargin lt width*4) then begin
    ;determine if there is enough blank space around the magnetogram to place a halo
    ;user has the choice to change the box size, or use another flux balancing method
       choice=0
       Print,'Insufficient blank space around magnetogram in box.'
       print,'Embedded magnetogram must fill no more than 75% of box'
       jump11: print,'Do you wish to'
       print,'(1) Choose another box size/scaling'
       print,'(2) Choose another flux balance method'
       read,choice
       if choice eq 1 then goto,jump4 else if choice eq 2 then goto,jump9 else begin
         print,'Invalid option! Please try again...'
         goto,jump11
       endelse
     endif  
     
     
     
     ;set up array containing the halo
     halo=fltarr(size,size)*0.
     halo(width/2:size-width/2-1,width/2:size-width/2-1)=1.
     halo(width*3/2:size-width*3/2-1,width*3/2:size-width*3/2-1)=0.
    
     n=n_elements(where(halo gt 0.5)) ;count number of pixels in halo
     
     for t=0,nt-1 do begin ;loop over each magnetogram
       dummy=mag_data(*,*,t)
       imb(t)=total(dummy) ;flux imbalance
       corr=imb(t)/float(n) ;pixel value in halo
       cor_mag(*,*,t) = dummy-corr*halo ;cor_mag is magnetogram + halo
     endfor
       
endif else begin
    print,'Invalid option! Please try again...'
    goto,jump9
endelse
      

      
      
fbskip: print,'--------------------------------------------------------------'
print,'Now showing a movie of the lower boundary condition:

;show movie of lower boundary time series
jumpmov:window,0,xs=size,ys=size
for t=0,nt-1 do begin
  tv,bytscl(cor_mag(*,*,t),-100,100)
  loadct,13,/silent
  if max([xstart,xstop,ystart,ystop]) gt 0.1 then begin
     plots,[xstart,xstart],[ystart,ystop],/device,linestyle=2,color=150
     plots,[xstop,xstop],[ystart,ystop],/device,linestyle=2,color=150
     plots,[xstart,xstop],[ystart,ystart],/device,linestyle=2,color=150
     plots,[xstart,xstop],[ystop,ystop],/device,linestyle=2,color=150
  endif
  
  if choice eq '2' then begin
     plots,[width/2,size-width/2],[width/2,width/2],/device,color=255,linestyle=2
     plots,[width/2,size-width/2.],[size-width/2,size-width/2],/device,color=255,linestyle=2
     plots,[width/2,width/2],[width/2,size-width/2],/device,color=255,linestyle=2
     plots,[size-width/2,size-width/2],[width/2,size-width/2],/device,color=255,linestyle=2
     
     plots,[width/2.*3.,size-width/2.*3.],[width/2.*3.,width/2.*3.],/device,color=255,linestyle=2
     plots,[width/2.*3.,size-width/2.*3.],[size-width/2.*3.,size-width/2.*3.],/device,color=255,linestyle=2
     plots,[width/2.*3.,width/2.*3.],[width/2.*3.,size-width/2.*3.],/device,color=255,linestyle=2
     plots,[size-width/2.*3.,size-width/2.*3.],[width/2.*3.,size-width/2.*3.],/device,color=255,linestyle=2
     
  endif
  loadct,0,/silent
  
  wait,0.05
endfor

answer='a'
jumptest:print,'Would you like to replay the movie? (y/n)'
read,answer
if answer eq 'y' then begin
  goto,jumpmov
endif else if answer ne 'n' then begin
  print, 'Invalid response! Please try again...'
  goto,jumptest
endif
print,''

wdelete,0
print,'--------------------------------------------------------------'

answer=0
;################# SELECT START & END FRAME ######################
print,'The magnetogram time series contains ',strtrim(nt,2),' frames.'
     Print,"Please refer to Window 5 for the active region's properties"
     print,"with time."
jump22: print,'Would you like to use'
print,'(1) The whole time series'
print,'(2) Only a portion of it'
read,answer
start=0
stop=0

;select start and end frames. Check that they are sensible choices
  if answer eq 2 then begin
     jump20: print,'Select the start frame'
     read,start
     if (start lt 0) or (start gt nt-2) then begin
        print,'Invalid start frame! Please choose again...'
        goto,jump20
     endif
     jump21: print,'Select the stop frame'
     read,stop
     if (stop lt 0) or (stop gt nt-1) then begin
        print,'Invalid stop frame! Please choose again...'
        goto,jump21
     endif
     if start gt stop then begin
        print,'Start frame cannot be after stop frame! Please choose again!'
        goto,jump20
     endif
     
     cor_mag=cor_mag(*,*,start:stop) ;remove discarded frames
     nt=stop-start+1 ;decrease nt
  endif else if answer ne 1 then begin
     print,'Invalid option! Please try again...'
     goto,jump22
  endif


end




