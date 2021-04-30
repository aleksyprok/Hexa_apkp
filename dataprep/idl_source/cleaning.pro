pro cleaning
@cleanup.blk
; This procedure takes the raw magnetogram (DATA) and produces a cleaned
; magnetogram from it (CLEANED). Various time-averaging options can be 
; chosen, along with low flux and small feature removal.



;################ TIME AVERAGING ############################
print,'--------------------------------------------------------------'
print,'Firstly time-averaging can be applied. Using this method each '
print,'cleaned frame is made from a linear combination of the raw'
print,'frames. This reduces noise and short-lived features in the '
print,'data.'
jump1: print,'Please select the time-averaging method to be applied:'
print,' (1) Boxcar'
print,' (2) Gaussian [Recommended]'
print,' (3) No time-averaging'
read,choice

if choice eq 3 then begin ; no smoothing

   cleaned=data

endif else if choice eq 1 then begin ; boxcar smoothing
   
   cleaned(*,*,0) = (data(*,*,0)+data(*,*,1))/2. ;first frame
   
   for t=1,nt-2 do begin
      cleaned(*,*,t) = (data(*,*,t-1)+data(*,*,t)+data(*,*,t+1))/3.;middle frames
   endfor
   
   cleaned(*,*,nt-1) = (data(*,*,nt-1)+data(*,*,nt-2))/2. ; end frame
   
endif else if choice eq 2 then begin ; Gaussian smoothing
   
   tau=2. ; smoothing time set to 2 frames
   
   for t=0,nt-1 do begin
       weight=exp(-(t-findgen(nt))^2/tau^2) ;construct weighting function
       gauss=0
       for k=0,nt-1 do begin
         gauss = gauss+ weight(k)*data(*,*,k) ; apply sum of weighted frames
       endfor
       
       cleaned(*,*,t)=gauss/total(weight) ;divide by sum of weights to get cleaned frame
   endfor
       
endif else begin
   print,'Invalid option... Please try again'
   goto, jump1
endelse


; ######################## LOW FIELD REMOVAL ###################
choice='a'
print,''
print,'--------------------------------------------------------------'
print,'Next low field removal can be applied. This sets all pixels'
print,'whose absolute values are less than 25 G to zero. An example'
print,'of this is displayed in Window 0 for frame '+strtrim(nt/2,2)+'.'
print,''

;preview low field removal
!p.background=255
window,0,xs=600,ys=300
frame=cleaned(*,*,nt/2)
!p.multi=[0,2,1]
pg_plotimage,bytscl(frame,-100,100),title='No low flux removal',charsize=1,/isotropic,color=0
zero=where(abs(frame) lt 25)
frame(zero)=0
pg_plotimage,bytscl(frame,-100,100),title='low flux removal',charsize=1,/isotropic,color=0

jump2: print,'Would you like to apply low field removal? (y[recommended]/n)'

read,choice

if choice eq 'y' then begin
   zero=where(abs(cleaned) lt 25) 
   cleaned(zero)=0
endif else if choice ne 'n' then begin
   print,'Invalid option... Please try again'
   goto, jump2
endif

wdelete,0
window,0,xs=600,ys=300




; ##################### ISOLATED FEATURE REMOVAL ####################
print,''
print,'--------------------------------------------------------------'
print,'The final cleaning operation that can be applied is isolated-'
print,'feature removal. This process removes small isolated features '
print,'in the magnetograms, such as cosmic ray hits. An example of'
print,'isolated-feature removal is displayed in Window 0 for frame'
print,strtrim(nt/2,2)+'.'

;preview isolated feature removal
frame=cleaned(*,*,nt/2)
pg_plotimage,bytscl(frame,-100,100),title='no isolated-feature removal',/isotropic,color=0
isolated,frame,nx,ny,0
pg_plotimage,bytscl(frame,-100,100),title='isolated-feature removal',/isotropic,color=0


!p.background=0
!p.multi=0

jump3: print,'Apply isolated feature removal? (y[recommended]/n)'
read,choice

if choice eq 'y' then begin
  isolated,cleaned,nx,ny,nt

endif else if choice ne 'n' then begin
   print,'Invalid option... Please try again'
   goto, jump3
endif

wdelete,0
print,''
print,''



end


