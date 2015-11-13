pro readdatap, file, x, y, z, v, time, step
;+
; read binary data
;-

print, file
NGH=2
openr, unit, file, /F77_UNFORMATTED, /GET_LUN 
time=0.d0 & step=0L
IMAX=0L & JMAX=0L & KMAX=0L & MMAX=0L
readu, unit, Time, Step
readu, unit, IMAX, JMAX, KMAX, MMAX
x=dblarr(IMAX+1+2*NGH)
y=dblarr(JMAX+1+2*NGH)
z=dblarr(KMAX+1+2*NGH)
v=dblarr(IMAX+1+2*NGH,JMAX+1+2*NGH,KMAX+1+2*NGH, MMAX+1)
readu, unit, X, Y, Z
readu, unit, V
FREE_LUN, unit    

x = x[NGH:IMAX+NGH]
y = y[NGH:JMAX+NGH]
z = z[NGH:KMAX+NGH]
v = v[NGH:IMAX+NGH,NGH:JMAX+NGH,NGH:KMAX+NGH,*]

end
