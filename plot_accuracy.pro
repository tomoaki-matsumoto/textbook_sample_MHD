;+
;
;-

; dir='DATA_green/'
;; dir='DATA_exp/'
dir='DATA_impl/'
bool_showExactsol=0
spawn, 'ls '+dir+'st*', files

!p.thick=2
!x.thick=2
!y.thick=2
!p.charsize=1
!x.charsize=1
!y.charsize=1
;; !p.font=0
!p.font=-1
colortable_linear, [[1, 240, 1, 1], [254, 0, 1, 1]],/hsv
set_color_index, 0, 255, 255, 255, /rgb
set_color_index, 255, 0, 0, 0, /rgb

error = dblarr(n_elements(files))
dtime = dblarr(n_elements(files))

for n = 0, n_elements(files)-1 do begin
   fn = files[n]
   NGH=1
   openr, unit, fn, /F77_UNFORMATTED, /GET_LUN 
   time=0.d0 & step=0L
   IMAX=0L & JMAX=0L & KMAX=0L & MMAX=0L
   readu, unit, Time, Step
   readu, unit, IMAX, JMAX, KMAX, MMAX
   x=dblarr(IMAX+1+2*NGH)
   y=dblarr(JMAX+1+2*NGH)
   z=dblarr(KMAX+1+2*NGH)
   u=dblarr(IMAX+1+2*NGH,JMAX+1+2*NGH,KMAX+1+2*NGH, MMAX+1)
   readu, unit, X, Y, Z
   readu, unit, U
   FREE_LUN, unit    

   t0=1.e-3
   eta = reform(u[*,1,1,1])
   B = reform(u[*,1,1,0])
   ex = 1/sqrt(4*!PI*eta*(Time+T0)) * exp(-x^2 / (4*eta*(Time+T0)) )
   error[n] = total(abs(B-ex))*(x(1)-x(0))
   dtime[n] = Time/Step
   print , Time
endfor
plot, dtime, error, /xlog, /ylog, xtitle=get_char('Delta')+'t', ytitle='L1 norm'
oplot, dtime, error, psym=4
oplot,[min(dtime),max(dtime)], (100*[min(dtime),max(dtime)])^2
oplot,[min(dtime),max(dtime)], (100*[min(dtime),max(dtime)])
end

