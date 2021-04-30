pro quick_analysis

dir='../hexaf90/run1' ;directory containing the hexa files

fname=dir+'/analysis.out'

openr,10,fname

t=[0]
h=[0.d0]
b2=[0]
bp2=[0]
ef=[0]
hdum=0.d0

while ~EOF(10) do begin
readf,10,tdum,hdum,b2dum,bp2dum,efdum
 t=[t,tdum]
 h=[h,hdum]
 b2=[b2,b2dum]
 bp2=[bp2,bp2dum]
 ef=[ef,efdum]
endwhile
close,10

n=n_elements(t)

t=t(1:n-1)
h=h(1:n-1)
b2=b2(1:n-1)
bp2=bp2(1:n-1)
ef=ef(1:n-1)

;remove,0,t,h,b2,bp2,ef
!p.multi=[0,1,3]
window,1,xs=800,ys=800
plot,t,h,charsize=2,xtitle='frame',ytitle='Helicity (MX !U2!N)'

;window,2
plot,t,b2,charsize=2,yrange=[0,max(b2)],xtitle='frame',ytitle='Magnetic energy (erg)'
oplot,t,bp2,linestyle=2

;window,3
plot,t,ef,charsize=2,xtitle='frame',ytitle='Free Magnetic energy (erg)'

!p.multi=0


stop
end