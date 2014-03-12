#include "config.h"
#ifndef FLUX_SCHEME_MHD
ERROR
#endif
module flux_eos
  use parameter
  implicit none
  private
  real(kind=DBL_KIND),parameter :: GAMMA = 5.d0/3.d0
!!$  real(kind=DBL_KIND),parameter :: GAMMA = 1.4d0
  real(kind=DBL_KIND),save :: Ch
  real(kind=DBL_KIND),parameter :: Cr = 0.18d0 ! Cr := cp^2/ch (for divB clean)
  public :: cyclecomp, flux, source_b, cflcond, v2u, u2v, GAMMA
contains
  !-----------------------------------------------------------------------
  ! rotates components of system equation
  !-----------------------------------------------------------------------
  function cyclecomp(ncrd, invert) result( mcycle )
    integer,intent(IN) :: ncrd
    integer,intent(IN),optional :: invert
    integer,dimension(MMIN:MMAX) :: mcycle
    integer,dimension(MX:MZ) :: mcycleV = (/ MVX, MVY, MVZ /)
    integer,dimension(MX:MZ) :: mcycleB = (/ MBX, MBY, MBZ /)
    integer :: m
    do m = MMIN, MMAX
       mcycle(m) = m
    enddo
    if ( present( invert ) ) then
       mcycle(MVX:MVX+size(mcycleV)-1) = cshift( mcycleV, -ncrd)
       mcycle(MBX:MBX+size(mcycleB)-1) = cshift( mcycleB, -ncrd)
    else
       mcycle(MVX:MVX+size(mcycleV)-1) = cshift( mcycleV, ncrd )
       mcycle(MBX:MBX+size(mcycleB)-1) = cshift( mcycleB, ncrd )
    endif
  end function cyclecomp
  !-----------------------------------------------------------------------
  ! get numerical flux in one dimension (Miyoshi & Kusano, 2005, JCP 208, 315-344)
  !-----------------------------------------------------------------------

  ! SWD(A, B) = 1.d0  if A >= B. otherwise 0.d0
#define SWD(A, B)  (0.5d0 + 0.5d0*sign(1.d0, (A)-(B)))
  ! SWR(A, B) = 1.d0 if A * B < 0, otherwise 0.d0
#define SWR(A, B)  (0.5d0 - 0.5d0*sign(1.d0, (A))*sign(1.d0, (B)))

  subroutine flux( ql, qr, f)
    use grid, only : get_cellWidth, Dtime
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: f   ! (OUT)
    real(kind=DBL_KIND) :: sw_deg_l, sw_deg_r, swl, swr, swal, swar, swml, swmr, swt
    real(kind=DBL_KIND),parameter :: eps = 1.D-6
    integer :: i,j,k, is,js,ks,ie,je,ke
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND) :: &
         rhol,ul,vl,wl,pl,bxl,byl,bzl,pbl,ptl,el, &
         rhoal,ual,val,wal,pal,byal,bzal,pbal,ptal,eal, &
         rhoml,uml,vml,wml,pml,byml,bzml,pbml,ptml,eml, &
         cf_l,csl2,cal2,caxl2,sl,sal,sql,psil,rhoil, &
         denom_l, &
                                ! 
         rhor,ur,vr,wr,pr,bxr,byr,bzr,pbr,ptr,er, &
         rhoar,uar,var,war,par,byar,bzar,pbar,ptar,ear, &
         rhomr,umr,vmr,wmr,pmr,bymr,bzmr,pbmr,ptmr,emr, &
         cf_r,csr2,car2,caxr2,sr,sar,sqr,psir,rhoir, &
         denom_r, &
                                ! 
         gm1i,sqrtpi4i,sqrtpi4,bxm,psim,sm,pta,sqlri,vm,wm,bym,bzm,gzm,signbxm,bxm2, &
         rho, u, v, w, bx, by, bz, p, e, pt

    h = get_cellWidth()
    Ch = (CFL) * minval(h) / Dtime / 3.d0

    is = max( lbound(ql,1), lbound(qr,1), lbound(f,1) )
    js = max( lbound(ql,2), lbound(qr,2), lbound(f,2) )
    ks = max( lbound(ql,3), lbound(qr,3), lbound(f,3) )
    ie = min( ubound(ql,1), ubound(qr,1), ubound(f,1) )
    je = min( ubound(ql,2), ubound(qr,2), ubound(f,2) )
    ke = min( ubound(ql,3), ubound(qr,3), ubound(f,3) )

    gm1i = 1.d0/(GAMMA - 1)
    sqrtpi4i = sqrt(PI4I)
    sqrtpi4 = sqrt(PI4)
    do k = ks, ke
       do j = js, je
          do i = is, ie
             ! -------------------------
             ! add-on for div B cleaning
             ! -------------------------
             bxl = ql(i,j,k,MBX)*sqrtpi4i
             bxr = qr(i,j,k,MBX)*sqrtpi4i
             psil = ql(i,j,k,MDB)*sqrtpi4i
             psir = qr(i,j,k,MDB)*sqrtpi4i
             bxm  = 0.5d0*(bxl+bxr - (psir-psil)/Ch)
             psim = 0.5d0*(psil+psir - (bxr-bxl)*Ch)
             f(i,j,k,MBX) = psim         * sqrtpi4
             f(i,j,k,MDB) = bxm*Ch**2    * sqrtpi4
             bxm2 = bxm**2
             ! --------------------
             ! U_L; left variables
             ! --------------------
             rhol = abs(ql(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
             ul = ql(i,j,k,MVX)
             vl = ql(i,j,k,MVY)
             wl = ql(i,j,k,MVZ)
             byl = ql(i,j,k,MBY)*sqrtpi4i
             bzl = ql(i,j,k,MBZ)*sqrtpi4i
             pl = abs(ql(i,j,k,MP))
             pbl = (bxm2 + byl**2 + bzl**2)*0.5d0
             ptl = pl + pbl
             el = rhol * (ul**2 + vl**2 + wl**2) * 0.5d0 + pbl + pl*gm1i
             ! ---------------------
             ! U_R; right variables
             ! ---------------------
             rhor = abs(qr(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
             ur = qr(i,j,k,MVX)
             vr = qr(i,j,k,MVY)
             wr = qr(i,j,k,MVZ)
             byr = qr(i,j,k,MBY)*sqrtpi4i
             bzr = qr(i,j,k,MBZ)*sqrtpi4i
             pr = abs(qr(i,j,k,MP))
             pbr = (bxm2 + byr**2 + bzr**2)*0.5d0
             ptr = pr + pbr
             er = rhor * (ur**2 + vr**2 + wr**2) * 0.5d0 + pbr + pr*gm1i
             ! -----------
             ! wave speed
             ! -----------
             rhoil = 1.d0/rhol
             rhoir = 1.d0/rhor
             csl2 = GAMMA * pl * rhoil ! sound speed
             csr2 = GAMMA * pr * rhoir ! sound speed
             cal2 = pbl * rhoil * 2.d0 ! alfven speed
             car2 = pbr * rhoir * 2.d0 ! alfven speed
             caxl2 = bxm2 * rhoil    ! alfven speed in x-direction
             caxr2 = bxm2 * rhoir    ! alfven speed in x-direction
             cf_l = sqrt( 0.5d0*( csl2 + cal2 + sqrt( (csl2 + cal2)**2 - 4.d0*csl2*caxl2) ) ) ! fast wave
             cf_r = sqrt( 0.5d0*( csr2 + car2 + sqrt( (csr2 + car2)**2 - 4.d0*csr2*caxr2) ) ) ! fast wave
             sl = min(ul-cf_l, ur-cf_r)
             sr = max(ul+cf_l, ur+cf_r)
             sm = ((sr - ur)*rhor*ur - (sl - ul)*rhol*ul - ptr + ptl)/((sr - ur)*rhor - (sl - ul)*rhol)
             rhoal = rhol * (sl - ul)/(sl - sm)
             rhoar = rhor * (sr - ur)/(sr - sm)
             sql = sqrt(rhoal)
             sqr = sqrt(rhoar)
             sal = sm - abs(bxm) / sql
             sar = sm + abs(bxm) / sqr
             ! ---------------------------
             ! switching flux
             ! ---------------------------
             swl  = SWD(sl, 0.d0)
             swal = SWR(sl, sal)
             swml = SWR(sal, sm)
             swmr = SWR(sm, sar)
             swar = SWR(sar, sr)
             swr  = SWD(0.d0, sr)
             swt = 1.d0/(swl+swr+swal+swar+swml+swmr)
             ! ---------------
             ! U*_L and U*_R
             ! ---------------
             pta = ((sr-ur)*rhor*ptl - (sl-ul)*rhol*ptr + rhol*rhor*(sr-ur)*(sl-ul)*(ur-ul))/((sr-ur)*rhor-(sl-ul)*rhol)

             denom_l = rhol*(sl-ul)*(sl-sm)-bxm2
             sw_deg_l = SWD(eps, abs(denom_l)) ! switch for degenerate
             denom_l = (1.d0-sw_deg_l) / (denom_l + sw_deg_l)
             ual = sm
             val = vl - byl*bxm*(sm-ul) * denom_l
             wal = wl - bzl*bxm*(sm-ul) * denom_l
             byal = byl*(rhol*(sl-ul)**2-bxm2) * denom_l
             bzal = bzl*(rhol*(sl-ul)**2-bxm2) * denom_l
             ptal = pta
             eal = ((sl-ul)*el-ptl*ul+ptal*sm+bxm*(ul*bxm+vl*byl+wl*bzl-ual*bxm-val*byal-wal*bzal))/(sl-sm)

             denom_r = rhor*(sr-ur)*(sr-sm)-bxm2
             sw_deg_r = SWD(eps, abs(denom_r))
             denom_r = (1.d0-sw_deg_r) / (denom_r + sw_deg_r)
             uar = sm
             var = vr - bxm*byr*(sm-ur) * denom_r
             war = wr - bxm*bzr*(sm-ur) * denom_r
             byar = byr*(rhor*(sr-ur)**2-bxm2) * denom_r
             bzar = bzr*(rhor*(sr-ur)**2-bxm2) * denom_r
             ptar = pta
             ear = ((sr-ur)*er-ptr*ur+ptar*sm+bxm*(ur*bxm+vr*byr+wr*bzr-uar*bxm-var*byar-war*bzar))/(sr-sm)
             ! --------------
             ! U**_L, U**_R
             ! --------------
             sqlri = 1.d0/(sql+sqr)
             signbxm = sign(1.d0,bxm)
             vm = (sql*val+sqr*var+(byar-byal)*signbxm)*sqlri
             wm = (sql*wal+sqr*war+(bzar-bzal)*signbxm)*sqlri
             bym = (sql*byar+sqr*byal+sql*sqr*(var-val)*signbxm)*sqlri
             bzm = (sql*bzar+sqr*bzal+sql*sqr*(war-wal)*signbxm)*sqlri

             rhoml = rhoal
             ptml = ptal
             uml = sm
             vml = vm
             wml = wm
             byml = bym
             bzml = bzm
             eml = eal - sql*(ual*bxm+val*byal+wal*bzal - uml*bxm-vml*byml-wml*bzml)*signbxm

             rhomr = rhoar
             ptmr = ptar
             umr = sm
             vmr = vm
             wmr = wm
             bymr = bym
             bzmr = bzm
             emr = ear + sqr*(uar*bxm+var*byar+war*bzar - umr*bxm-vmr*bymr-wmr*bzmr)*signbxm
             ! -------------
             ! switching U
             ! -------------
             rho = (rhol*swl + rhor*swr + rhoal*swal + rhoar*swar + rhoml*swml + rhomr*swmr)*swt
             u = (ul*swl + ur*swr + ual*swal + uar*swar + uml*swml + umr*swmr)*swt
             v = (vl*swl + vr*swr + val*swal + var*swar + vml*swml + vmr*swmr)*swt
             w = (wl*swl + wr*swr + wal*swal + war*swar + wml*swml + wmr*swmr)*swt
             bx = bxm
             by = (byl*swl + byr*swr + byal*swal + byar*swar + byml*swml + bymr*swmr)*swt
             bz = (bzl*swl + bzr*swr + bzal*swal + bzar*swar + bzml*swml + bzmr*swmr)*swt
             pt = (ptl*swl + ptr*swr + ptal*swal + ptar*swar + ptml*swml + ptmr*swmr)*swt
             e = (el*swl + er*swr + eal*swal + ear*swar + eml*swml + emr*swmr)*swt
             ! --------
             ! flux
             ! --------
             f(i,j,k,MRHO) = rho*u
             f(i,j,k, MVX) = rho*u**2 + pt - bx**2
             f(i,j,k, MVY) = rho*v*u - bx*by
             f(i,j,k, MVZ) = rho*w*u - bx*bz
             f(i,j,k, MBY) = (by*u - bx*v)*sqrtpi4
             f(i,j,k, MBZ) = (bz*u - bx*w)*sqrtpi4
             f(i,j,k,  MP) = (e + pt)*u - bx*(u*bx+v*by+w*bz)
          enddo
       enddo
    enddo
  end subroutine flux
  !-----------------------------------------------------------------------
  ! source term due to div B
  !-----------------------------------------------------------------------
  subroutine source_b(f, w, dt)
    use util
    use grid, only : get_ds, V
    real(kind=DBL_KIND),dimension(IMINGH:,JMINGH:,KMINGH:,MMIN:,MX:) :: f  ! (IN)
    real(kind=DBL_KIND),dimension(IMINGH:,JMINGH:,KMINGH:,MMIN:) :: w  ! (OUT)
    real(kind=DBL_KIND),intent(IN) :: dt
    real(kind=DBL_KIND),dimension(MX:MZ) :: ds
    integer :: i,j,k,n, mb, io,jo,ko
    real(kind=DBL_KIND) :: divbdv, dtchcr
    integer,dimension(MMIN:MMAX) :: mcycle

    ds = get_ds()
    dtchcr = dt*Ch/Cr
    ! f(i,j,k,MBX,MX) = Psi @ x-surface
    ! f(i,j,k,MBY,MY) = Psi @ y-surface
    ! f(i,j,k,MBZ,MZ) = Psi @ z-surface
    ! f(i,j,k,MDB,MX) = ch^2 Bxm

    do n = MX, MX+NDIM-1
       mcycle = cyclecomp(n)
       mb = mcycle(MBX)
       call util_arroffset(n, io,jo,ko)
       do k=Kmin,Kmax
          do j=Jmin,Jmax
             do i=Imin,Imax
                divbdv = (f(i,j,k,MDB,n) - f(i-io,j-jo,k-ko,MDB,n))*ds(n)/(PI4*Ch**2)
                w(i,j,k,MVX)=w(i,j,k,MVX)-dt*divbdv*V(i,j,k,MBX)
                w(i,j,k,MVY)=w(i,j,k,MVY)-dt*divbdv*V(i,j,k,MBY)
                w(i,j,k,MVZ)=w(i,j,k,MVZ)-dt*divbdv*V(i,j,k,MBZ)
                w(i,j,k,MP)=w(i,j,k,MP) -dt*V(i,j,k,mb)*(f(i,j,k,mb,n)-f(i-io,j-jo,k-ko,mb,n))*ds(n)*PI4I
             enddo
          enddo
       enddo
    enddo
    do k=Kmin,Kmax
       do j=Jmin,Jmax
          do i=Imin,Imax
             w(i,j,k,MDB)=w(i,j,k,MDB)*exp(-dtchcr)
          end do
       end do
    end do
  end subroutine source_b
  !-----------------------------------------------------------------------
  ! find dt according to CFL condtion
  !-----------------------------------------------------------------------
  subroutine cflcond
    use grid
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: vx, vy, vz, rho, p, bx, by, bz
    real(kind=DBL_KIND),dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: ca
    h = get_cellWidth()
    rho => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MRHO)
    p   => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MP)
    vx  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MVX)
    vy  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MVY)
    vz  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MVZ)
    bx  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MBX)
    by  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MBY)
    bz  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MBZ)
    ca= sqrt( ( GAMMA*p +(bx**2+by**2+bz**2)*PI4I )/rho ) ! fast wave
#if defined(DIRECTIONAL_SPLIT)
    if (NDIM == 3) then
       Dtime = (CFL) / max( &
            maxval(abs(vx)+ca)/h(MX), &
            maxval(abs(vy)+ca)/h(MY), &
            maxval(abs(vz)+ca)/h(MZ) )
    elseif (NDIM == 2) then
       Dtime = (CFL) / max(  &
            maxval(abs(vx)+ca)/h(MX), &
            maxval(abs(vy)+ca)/h(MY) )
    elseif (NDIM == 1) then
       Dtime = (CFL)*minval( h(MX)/(abs(vx)+ca) )
    else
       print *, '*** error'
    endif
#elif defined(DIRECTIONAL_UNSPLIT)
    if (NDIM == 3) then
       Dtime = (CFL) / maxval( &
            (abs(vx)+ca)/h(MX) + &
            (abs(vy)+ca)/h(MY) + &
            (abs(vz)+ca)/h(MZ) )
    elseif (NDIM == 2) then
       Dtime = (CFL) / maxval(  &
            (abs(vx)+ca)/h(MX) + &
            (abs(vy)+ca)/h(MY) )
    elseif (NDIM == 1) then
       Dtime = (CFL)*minval( h(MX)/(abs(vx)+ca) )
    else
       print *, '*** error'
    endif
#else
    ERROR
#endif
  end subroutine cflcond
  !-----------------------------------------------------------------------
  ! convert primitive variables u to conservative variable u
  !-----------------------------------------------------------------------
  subroutine v2u(q,u)
    use grid, only : get_dv
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: q !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: u !(OUT)
    real(kind=DBL_KIND) :: dv, pi8i
    dv = get_dv()
    pi8i = PI4I/2.d0
    u(:,:,:,MRHO)=q(:,:,:,MRHO)*dv
    u(:,:,:,MVX)=q(:,:,:,MVX)*u(:,:,:,MRHO)
    u(:,:,:,MVY)=q(:,:,:,MVY)*u(:,:,:,MRHO)
    u(:,:,:,MVZ)=q(:,:,:,MVZ)*u(:,:,:,MRHO)
    u(:,:,:,MBX)= q(:,:,:,MBX)*dv
    u(:,:,:,MBY)= q(:,:,:,MBY)*dv
    u(:,:,:,MBZ)= q(:,:,:,MBZ)*dv
    u(:,:,:,MP) =dv*(q(:,:,:,MP)/(GAMMA-1.d0) &
         +q(:,:,:,MRHO)*(q(:,:,:,MVX)**2+q(:,:,:,MVY)**2+q(:,:,:,MVZ)**2)/2.d0 &
         +(q(:,:,:,MBX)**2+q(:,:,:,MBY)**2+q(:,:,:,MBZ)**2)*pi8i)
    u(:,:,:,MDB) = q(:,:,:,MDB)*dv
  end subroutine v2u
  !-----------------------------------------------------------------------
  ! convert conservative variables u to primitive variable u
  !-----------------------------------------------------------------------
  subroutine u2v(u,q)
    use grid, only : get_dv
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: q !(OUT)
    real(kind=DBL_KIND) :: dv, pi8i
    dv = get_dv()
    pi8i = PI4I/2.d0
    q(:,:,:,MRHO) =u(:,:,:,MRHO)/dv
    q(:,:,:,MVX) =u(:,:,:,MVX)/u(:,:,:,MRHO)
    q(:,:,:,MVY) =u(:,:,:,MVY)/u(:,:,:,MRHO)
    q(:,:,:,MVZ) =u(:,:,:,MVZ)/u(:,:,:,MRHO)
    q(:,:,:,MBX) =u(:,:,:,MBX)/dv
    q(:,:,:,MBY) =u(:,:,:,MBY)/dv
    q(:,:,:,MBZ) =u(:,:,:,MBZ)/dv
    q(:,:,:,MP) = ( u(:,:,:,MP)/dv -0.5d0*q(:,:,:,MRHO) &
         *(q(:,:,:,MVX)**2+q(:,:,:,MVY)**2+q(:,:,:,MVZ)**2) &
         -(q(:,:,:,MBX)**2+q(:,:,:,MBY)**2+q(:,:,:,MBZ)**2)*pi8i &
         )*(GAMMA-1.0d0)
    q(:,:,:,MDB) = u(:,:,:,MDB)/dv
  end subroutine u2v
  !-----------------------------------------------------------------------
  ! return sound speed (c) when density (rho) is given
  !-----------------------------------------------------------------------
  subroutine get_cs(c, rho, p)
    real(kind=DBL_KIND),intent(OUT) :: c
    real(kind=DBL_KIND),intent(IN) :: rho, p
    c = sqrt(abs(GAMMA*p/rho))
  end subroutine get_cs


end module flux_eos
