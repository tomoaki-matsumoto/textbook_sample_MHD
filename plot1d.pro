;+
;
;-

dir='DATA'
;; dir = 'DATA_unsplit_predictor'
spawn, 'ls '+dir, r
fn = r[(size(r))[1]-1]
;; fn = r[0]
;; fn = 'st000150.d'

file = dir +'/'+fn
readdatap, file, x, y, z, v, time, step

; sound speed
gamma = 5.d0/3.d0
cs = sqrt(gamma*v[*,*,*,4]/v[*,*,*,0])
v[*,*,*,1]=v[*,*,*,1]/cs

if (size(v))[1] eq n_elements(x) then begin
   mplot=1
   plot, x, v[*,0,0,mplot],yrange=[min(v[*,0,0,mplot]), max(v[*,0,0,mplot])] ; , title='1st Order'
   oplot, x, v[*,0,0,mplot],psym=4
;   oplot,[1,1]*time, [0,3]
endif else if size(v,/n_dim) eq 2 then begin
   mplot=0
   tvcn, v[*,*,0,mplot], v[*,*,0,mplot], x, y, /noby,/asp, level=[0]
   ;; temp = v[*,*,0,7]/v[*,*,0,0]
   ;; tvcn, temp, temp, x, y, /noby,/asp, level=[0]
   ;; beta = (v[*,*,0,4]^2+v[*,*,0,5]^2+v[*,*,0,6]^2)/v[*,*,0,7]
   ;; tvcn, beta, beta, x, y, /noby,/asp, level=[0]
endif 

end
