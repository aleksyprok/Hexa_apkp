pro isolated ,cleaned,nx,ny,nt
;Removes isolated pixels from a magnetogram

for t=0,nt-1 do begin
    ;clean=cleaned(*,*,t)
    for i=1,nx-2 do begin
      for j=1,ny-2 do begin
	surround=cleaned(i-1:i+1,j-1:j+1,t) ;group of 8 neighbouring pixels
	surround=surround*surround(1,1)
	surround(1,1)=0
	test=where(surround gt 0) ; where neighbours have the same sign as the centre pixel
	if n_elements(test) lt 4 then begin
	   cleaned(i,j,t)=0 ;set centre pixel with too few neighbours of same sign to zero
	endif
      endfor
    endfor
endfor


;set pixels around the edge of the magnetogram to zero
cleaned(0,*,*)=0
cleaned(nx-1,*,*)=0
cleaned(*,0,*)=0
cleaned(*,ny-1,*)=0


end