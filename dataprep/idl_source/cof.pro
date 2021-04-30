pro cof, data,nx,ny,x,y
;calculates centre of flux (x,y) of a magnetogram (DATA) of dimenions (nx,ny)
x=0.
y=0.

for j=0,ny-1 do begin
  for i=0,nx-1 do begin
     x=x+i*data(i,j)
     y=y+j*data(i,j)
  endfor
endfor
x=x/total(data)
y=y/total(data)

end
    