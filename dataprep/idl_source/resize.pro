pro resize, res,embed=embed,resize=resize
@cleanup.blk
;Takes the magnetogram time series and inserts it into the computational box.
;The magnetograms can either be scaled to fll the box, be placed inside the
;box as-is, or have an arbitrary scaling inside the box.


mag_data=fltarr(res,res,nt)

wdelete,0

scale=float(res)/float(max([nx,ny])) ;determine how much the magnetogram must be scaled by for it to fill box
  nxnew=fix(nx*scale) ;scaled nx
  nynew=fix(ny*scale) ;scaled ny
if keyword_set(resize) then begin  ;Apply scaling to fill the box
  dummy=congrid(cleaned,nxnew,nynew,nt,/interp)
endif else if keyword_set(embed) then begin ;Apply arbitrary scaling
  choice=0
  jump12: print,'Would you like to:'
  print,'(1) Insert the magnetograms into the centre of the box at'
  print,'    their current size?'
  print,'(2) Custom scaling'
  read, choice
  if choice eq 1 then begin ; no scaling
     if max([nx,ny]) gt size then begin
        print,'Magnetogram larger than box. Please choose another option!'
        goto,jump12
     endif
     dummy=cleaned
     nxnew=nx
     nynew=ny     
  endif else if choice eq 2 then begin ;arbitrary scaling
     scale=5.
     jump7: print,'Enter size of magnetogram relative to size of box (range=[0,1])'
     read,scale
     if scale gt 1. then begin
        print, 'Magnetogram cannot be larger than box!'
        print, 'Please try again...'
        goto,jump7
     endif else if scale le 0 then begin
        print, 'Magnetogram cannot have negative/zero size!'
        print, 'Please try again...'
        goto,jump7
     endif else begin
        nxnew=fix(nxnew*scale)
        nynew=fix(nynew*scale)
        dummy=congrid(cleaned,nxnew,nynew,nt,/interp)
     endelse
  endif else begin
     print,'Invalid option! Please try again...'
     goto,jump12
  endelse
endif else begin
  print,'ERROR! Invalid keyword'
  stop
endelse


l0=float(nx)*deltx*float(res)/float(nxnew)/6. ;length of one unit in Hexa


;embed dummy (array containing scaled magnetograms) into the computational box

xmargin=res-nxnew
ymargin=res-nynew

if min([xmargin,ymargin]) lt 0 then begin ;check that resized magnetogram isn't bigger than box
  print,'data too big'
  stop
endif

;determine x coordinates of the edges of the scaled magnetogram within the box 
if even(xmargin) then begin
   xstart=xmargin/2
   xstop=res-xmargin/2-1
endif else begin
   xstart=xmargin/2
   xstop=res-xmargin/2-2
endelse
;determine y coordinates of the edges of the scaled magnetogram within the box 
if even(ymargin) then begin
   ystart=ymargin/2
   ystop=res-ymargin/2-1
endif else begin
   ystart=ymargin/2
   ystop=res-ymargin/2-2
endelse

;place scaled magnetograms in box
mag_data(xstart:xstop,ystart:ystop,*)=dummy
window,0,xs=res,ys=res
tv,bytscl(mag_data(*,*,nt/2),-100,100)
loadct,13,/silent
plots,[xstart,xstart],[ystart,ystop],/device,linestyle=2,color=150
plots,[xstop,xstop],[ystart,ystop],/device,linestyle=2,color=150
plots,[xstart,xstop],[ystart,ystart],/device,linestyle=2,color=150
plots,[xstart,xstop],[ystop,ystop],/device,linestyle=2,color=150
loadct,0,/silent
end



