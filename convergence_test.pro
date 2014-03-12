;+
;
;-

logfile = 'MUSCL3RK3.log'
;; logfile = 'MUSCL3RK3splt.log'
;; logfile = 'MUSCL3K3.log'
;; logfile = 'MUSCL2RK2.log'
;; logfile = 'MUSCL2PC2.log'
;; logfile = 'NoneEuler.log'
;; logfile = 'MUSCL3PC2.log'
;; logfile = 'AdvMUSCL3RK3.log'
;; logfile = 'AdvMUSCL2RK2.log'
;; logfile = 'AdvMUSCL3RK2.log'
;; logfile = 'AdvNoneEuler.log'

;; logfile = 'test.log'

logdir = 'convergence_test'

logfile = logdir +'/'+ logfile
nline = wc(logfile)
nz = lonarr(nline)
l1norm = dblarr(nline)
l2norm = dblarr(nline)
linfnorm = dblarr(nline)
t_nz=0L
t_norm=dblarr(3)
openr, unit, logfile, /GET_LUN
n = 0L
while ~ eof(unit) do begin
   readf, unit, t_nz, t_norm, format='(I5,3E12.5)'
   nz[n] = t_nz
   l1norm[n] = t_norm[0]
   l2norm[n] = t_norm[1]
   linfnorm[n] = t_norm[2]
   n = n + 1
endwhile
free_lun, unit
yrange=[min(l1norm)<min(l2norm)<min(linfnorm), max(linfnorm)>max(l2norm)>max(linfnorm)]
plot, nz, l1norm, yrange=yrange, /xlog, /ylog
oplot, nz, l1norm, psym=4
oplot, nz, l2norm, linestyle=1
oplot, nz, l2norm, psym=4
oplot, nz, linfnorm, linestyle=2
oplot, nz, linfnorm, psym=4

np=n_elements(nz)-2
norm=l1norm
; np=0
oplot, nz, norm[np]*(float(nz[np])/nz)
oplot, nz, norm[np]*(float(nz[np])/nz)^2
oplot, nz, norm[np]*(float(nz[np])/nz)^3

end

