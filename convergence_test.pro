;+
;
;-

;;;; Hydrodynamical wave
;; logfile = 'NoneEuler.log'
;; logfile = 'MUSCL2PC2.log'
;; logfile = 'MUSCL2RK2.log'
;; logfile = 'MUSCL3PC2.log'
;; logfile = 'MUSCL3RK3.log'
;; logfile = 'MUSCL3K3.log'

;;;; scalar advection
;; logfile = 'AdvNoneEuler.log'
;; logfile = 'AdvMUSCL2Euler.log'
;; logfile = 'AdvMUSCL2PC2.log'
logfile = 'AdvSuperBee2PC2.log'
;; logfile = 'AdvMUSCL2RK2.log'
;; logfile = 'AdvMUSCL3PC2.log'
;; logfile = 'AdvMUSCL3RK3.log'
;; logfile = 'AdvMUSCL3K3.log'

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
yrange=yrange * [0.1^0.1, 10^0.1]
xrange = [min(nz),max(nz)] * [0.1^0.1, 10^0.1]
plot, nz, l1norm, psym=usrsymbol(-10), xrange=xrange, yrange=yrange, /xlog, /ylog, xtitle='Number of cells', ytitle='Norms of error', xstyle=1, ystyle=1
oplot, nz, l2norm, psym=-7
oplot, nz, linfnorm, psym=-4


np=n_elements(nz)-1
;; norm=l1norm
normmin = l2norm[np]<l1norm[np]<linfnorm[np]
normmax = l2norm[np]>l1norm[np]>linfnorm[np]
; np=0
oplot, nz, 10^0.1* normmax*(float(nz[np])/nz),  linestyle=2
oplot, nz, 0.1^0.2* normmin*(float(nz[np])/nz)^2,linestyle=2
;; oplot, nz, 0.1^0.1* normmin*(float(nz[np])/nz)^3,linestyle=2

;; oplot, nz, normmax*(float(nz[np])/nz)^1.4,linestyle=2

lines=[0, 0, 0, 2]
psym = [-4, -7, -10, 0]
items = ['L'+get_char('infty'), 'L2', 'L1', 'O(n),n=1,2']
astron_legend,items,linestyle=lines, psym=psym, /right, /top


end

