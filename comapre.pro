dirs = [$
       'DATA_exp1st/', $
       'DATA_exp2nd/', $
       'DATA_fullimpl/', $
       'DATA_crankNicolson/']

colortable_linear, [[1, 240, 1, 1], [254, 0, 1, 1]],/hsv
set_color_index, 0, 255, 255, 255, /rgb
set_color_index, 255, 0, 0, 0, /rgb

for n = 0L, n_elements(dirs)-1 do begin
   dir = dirs[n]
   spawn, 'ls '+dir+'st*', files
   fn = files[3]
   NGH=1
   openr, unit, fn, /F77_UNFORMATTED, /GET_LUN 
   time=0.d0 & step=0L
   IMAX=0L & JMAX=0L & KMAX=0L & MMAX=0L
   readu, unit, Time, Step
   readu, unit, IMAX, JMAX, KMAX, MMAX
   if n eq 0 then begin
      x=dblarr(IMAX+1+2*NGH, n_elements(dirs))
      y=dblarr(JMAX+1+2*NGH, n_elements(dirs))
      z=dblarr(KMAX+1+2*NGH, n_elements(dirs))
      u=dblarr(IMAX+1+2*NGH,JMAX+1+2*NGH,KMAX+1+2*NGH, MMAX+1, n_elements(dirs))
      xx=dblarr(IMAX+1+2*NGH)
      yy=dblarr(JMAX+1+2*NGH)
      zz=dblarr(KMAX+1+2*NGH)
      uu=dblarr(IMAX+1+2*NGH,JMAX+1+2*NGH,KMAX+1+2*NGH, MMAX+1)
   endif
   readu, unit, xx, yy, zz
   readu, unit, uu
   x[*,n] = xx
   y[*,n] = yy
   z[*,n] = zz
   u[*,*,*,*,n] = uu
   print, dir, time
   FREE_LUN, unit    
endfor

plot, [min(x),max(x)],[min(u[*,*,*,0,*]),max(u[*,*,*,0,*])],/nodata, xtitle='x', ytitle='B', xstyle=1, $
;;       /ylog
      xrange=1+[-1,1]*0.1,yrange=.9+[-1,1]*0.1

for n = 0L, n_elements(dirs)-1 do begin
   color = fix(253*float(n)/(n_elements(dirs)-1))+1
   B = reform(u[*,1,1,0,n])
;;    B = abs(reform(u[*,1,1,0,n]) - reform(u[*,1,1,0,1]))
   oplot,x,B, color=color
   xyouts, 0.8, 0.9 - float(n)/(n_elements(dirs)-1) * 0.2, strcompress(dirs[n]), color=color,/norm
endfor

end
