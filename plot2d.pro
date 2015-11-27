;+
;
;-
dir='DATA'

spawn, 'ls '+dir, r
;; fn = r[0]
fn = r[(size(r))[1]-1]
file = dir +'/'+fn
readdatap, file, x, y, z, v, time, step
v=reform(v[*,*,0,0])

level=cont_level(v, nlevel=30)
cn, v, x, y,/asp, xtitle='x', ytitle='y',level=level


end

