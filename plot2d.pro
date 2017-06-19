;+
;
;-
dir='DATA'

spawn, 'ls '+dir, r
;; fn = r[0]
fn = r[(size(r))[1]-1]
file = dir +'/'+fn
readdatap, file, x, y, z, v, time, step

;; v=reform(v[*,*,0,0])
;; level=cont_level(v, nlevel=30)
;; cn, rho, x, y,/asp, xtitle='x', ytitle='y',level=level


rho = v[*,*,0,0]
vx = v[*,*,0,1]
vy = v[*,*,0,2]
vv = sqrt(vx^2+vy^2)
tvcnve,rho,vv,vx,vy,x, y,/asp, xtitle='x', ytitle='y',/noby,level=[1,2,3]


end

