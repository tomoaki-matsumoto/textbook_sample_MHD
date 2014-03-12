;+
;
;-


dir='DATA'
;; dir = 'DATA_unsplit_predictor'
spawn, 'ls '+dir, r
fn = r[(size(r))[1]-1]
; fn = r[0]
;; fn = 'st000150.d'

fn = dir +'/'+ fn
print, fn
NGH=2
openr, unit, fn, /F77_UNFORMATTED, /GET_LUN 
time=0.d0 & step=0L
IMAX=0L & JMAX=0L & KMAX=0L & MMAX=0L
readu, unit, Time, Step
readu, unit, IMAX, JMAX, KMAX, MMAX
x=dblarr(IMAX+1+2*NGH)
y=dblarr(JMAX+1+2*NGH)
z=dblarr(KMAX+1+2*NGH)
q=dblarr(IMAX+1+2*NGH,JMAX+1+2*NGH,KMAX+1+2*NGH, MMAX+1)
readu, unit, X, Y, Z
readu, unit, Q
FREE_LUN, unit    

x = x[NGH:IMAX+NGH]
y = y[NGH:JMAX+NGH]
z = z[NGH:KMAX+NGH]
q = q[NGH:IMAX+NGH,NGH:JMAX+NGH,NGH:KMAX+NGH,*]

if n_elements(q) eq (IMAX+1)*(MMAX+1) then begin
   mplot=0
   plot, x, q[*,0,0,mplot],yrange=[min(q[*,0,0,mplot]), max(q[*,0,0,mplot])] ; , title='1st Order'
   oplot, x, q[*,0,0,mplot],psym=4
   oplot,[1,1]*time, [0,3]
endif else if size(q,/n_dim) eq 2 then begin
   mplot=0
   tvcn, q[*,*,0,mplot], q[*,*,0,mplot], x, y, /noby,/asp, level=[0]
   ;; temp = q[*,*,0,7]/q[*,*,0,0]
   ;; tvcn, temp, temp, x, y, /noby,/asp, level=[0]
   ;; beta = (q[*,*,0,4]^2+q[*,*,0,5]^2+q[*,*,0,6]^2)/q[*,*,0,7]
   ;; tvcn, beta, beta, x, y, /noby,/asp, level=[0]
endif 

end
