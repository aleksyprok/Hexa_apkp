

t=[0]
dt=[0]
e=[0]
sz=[0]
q=[0]
eta=[0]
eta4=[0]
j=[0.d0]
jdum=0.d0
openr,10,'diagnostic'

while ~EOF(10) do begin
readf,10,tdum,dtdum,edum,sdum,qdum,etadum,eta4dum,jdum
t=[t,tdum]
dt=[dt,dtdum]
e=[e,edum]
sz=[sz,sdum]
q=[q,qdum]
eta=[eta,etadum]
eta4=[eta4,eta4dum]
j=[j,jdum]
endwhile

close,10



base=total(sz*dt,/cumulative)
diss=total(q*dt,/cumulative)
ohmic=total(eta*dt,/cumulative)
hyper=total(eta4*dt,/cumulative)


diss2=diss+hyper+ohmic





window,3,xs=500,ys=1000
!p.multi=[0,1,4]

plot,t,dt,xtitle='time',ytitle='timestep (s)',psym=3,charsize=2
oplot,[0,max(t)],[60,60],linestyle=2



plot,t,e,yrange=[0,max(e+diss2)],xtitle='time',ytitle='Energies (erg)',charsize=2
oplot,t,e+diss2,linestyle=2


plot,t,q+eta+eta4,xtitle='time',ytitle='Dissipation',charsize=2
oplot,t,q,linestyle=2
oplot,t,eta,linestyle=3
oplot,t,eta4,linestyle=4

plot,t,j,yrange=[0,max(j)],xtitle='time',ytitle='volume integrated current',charsize=2



!p.multi=0


print,'--------Summary---------'
print,'Mean timestep: '+strtrim(mean(dt),2)+' s'
print,'Mean dissipation: '+strtrim(mean(q),2)+' erg/s'
print,'Total energy dissipated: '+strtrim(total(q*dt),2)+' erg'



 stop
 
end