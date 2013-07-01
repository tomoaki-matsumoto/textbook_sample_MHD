#include "config.h"
module flux_eos
  use parameter
  implicit none
  private
  real(kind=DBL_KIND),parameter :: GAMMA = 5.d0/3.d0
!!$  real(kind=DBL_KIND),parameter :: GAMMA = 1.4d0
  public :: cyclecomp, flux, cflcond, v2u, u2v
contains
  !-----------------------------------------------------------------------
  ! rotates components of system equation
  !-----------------------------------------------------------------------
  function cyclecomp(ncrd, invert) result( mcycle )
    integer,intent(IN) :: ncrd
    integer,intent(IN),optional :: invert
    integer,dimension(Mmin:Mmax) :: mcycle
    integer,dimension(MX:MZ) :: mcycle3 = (/ MVX, MVY, MVZ /)
    integer :: m
    do m = Mmin, Mmax
       mcycle(m) = m
    enddo
    if ( present( invert ) ) then
       mcycle(MVX:MVX+size(mcycle3)-1) = cshift( mcycle3, -ncrd)
    else
       mcycle(MVX:MVX+size(mcycle3)-1) = cshift( mcycle3, ncrd )
    endif
  end function cyclecomp
  !-----------------------------------------------------------------------
  ! get numerical flux in one dimension (RoeM2; Kim et al. 2003, JCP, 185, 342)
  !-----------------------------------------------------------------------
  ! macro for entropy condition
#define ELMOD(EL,ELBAR,EL1, EL2) \
  eps=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-eps) ;\
  x2=1.d0-x1 ;\
  EL= x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(eps+x1)+eps)

  subroutine flux( ql, qr, f)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: f   ! (OUT)
    integer :: i,j,k, is,js,ks,ie,je,ke
    real(kind=DBL_KIND) :: &
         rhol,ul,vl,wl,pl,hl,el,cl, &
         rhor,ur,vr,wr,pr,hr,er,cr, &
         drho, drhou, drhov, drhow, drhoh, dp, dh, du, dv, dw, &
         sql,sqr,sqa,rhob,ub,vb,wb,hb,qb2,cb2,cb,ub2,vb2,wb2, &
         gm1,mn, b1, b2, b1b2, db12, pratio, ff, gg, rbdq, &
         fl1,fl2,fl3,fl4,fl5, &
         fr1,fr2,fr3,fr4,fr5 
!!$    real(kind=DBL_KIND) :: el4c, a7c, a8c, a9c, w1c 
    real(kind=DBL_KIND) :: eps, x1, x2
    is = max( lbound(ql,1), lbound(qr,1), lbound(f,1) )
    js = max( lbound(ql,2), lbound(qr,2), lbound(f,2) )
    ks = max( lbound(ql,3), lbound(qr,3), lbound(f,3) )
    ie = min( ubound(ql,1), ubound(qr,1), ubound(f,1) )
    je = min( ubound(ql,2), ubound(qr,2), ubound(f,2) )
    ke = min( ubound(ql,3), ubound(qr,3), ubound(f,3) )

    do k = ks, ke
       do j = js, je
          do i = is, ie
             gm1 = GAMMA - 1
             ! --------------
             ! left variables
             ! --------------
             rhol = abs(ql(i,j,k,MRHO)) ! 袖の未定義値参照でお亡くらないように。
             ul = ql(i,j,k,MVX)
             vl = ql(i,j,k,MVY)
             wl = ql(i,j,k,MVZ)
             pl = abs(ql(i,j,k,MP))
             el = rhol * (ul**2 + vl**2 + wl**2) / 2 + pl / gm1
             hl = ( el + pl )/rhol
             cl = sqrt( GAMMA * pl / rhol )
             ! --------------
             ! right variables
             ! --------------
             rhor = abs(qr(i,j,k,MRHO))
             ur = qr(i,j,k,MVX)
             vr = qr(i,j,k,MVY)
             wr = qr(i,j,k,MVZ)
             pr = abs(qr(i,j,k,MP))
             er = rhor * (ur**2 + vr**2 + wr**2) / 2 + pr / gm1
             hr = ( er + pr )/rhor
             cr = sqrt( GAMMA * pr / rhor )
             ! -------
             ! delta Q
             ! -------
             drho = rhor - rhol
             drhou = rhor * ur - rhol * ul
             drhov = rhor * vr - rhol * vl
             drhow = rhor * wr - rhol * wl
             drhoh = rhor * hr - rhol * hl
             du = ur - ul
             dv = vr - vl
             dw = wr - wl
             dp = pr - pl
             dh = hr - hl
             ! ----------------------
             ! intermidiate variables
             ! ----------------------
             sql = sqrt(rhol)
             sqr = sqrt(rhor)
             sqa = sql + sqr
             rhob = sql * sqr
             ub = (sql * ul + sqr * ur) / sqa
             vb = (sql * vl + sqr * vr) / sqa
             wb = (sql * wl + sqr * wr) / sqa
             hb = (sql * hl + sqr * hr) / sqa
             ub2 = ub**2
             vb2 = vb**2
             wb2 = wb**2
             qb2 = ub2 + vb2 + wb2
             cb2 = gm1*(hb - qb2/2)
             cb  =  sqrt( cb2 )
             mn = abs(ub/cb)         !Mach number

             ! 
             b1 = max(0.d0, ub+cb, ur+cb)
             b2 = min(0.d0, ub-cb, ul-cb)
             b1b2 = b1 * b2
             db12 = b1 - b2
             pratio = min(pr/pl, pl/pr)
             if (mn == 0.d0) then
                ff = 1.d0
             else
                ff = mn**(1.d0-pratio)
             endif
             gg = ff
             rbdq = drho - ff * dp / cb2
             ! ---------
             ! FL and FR
             ! ---------
             fl1 = rhol * ul
             fl2 = fl1 * ul + pl
             fl3 = fl1 * vl
             fl4 = fl1 * wl
             fl5 = ul * ( el + pl )

             fr1 = rhor * ur
             fr2 = fr1 * ur + pr
             fr3 = fr1 * vr
             fr4 = fr1 * wr
             fr5 = ur * ( er + pr )
             ! -----------------
             ! intermidiate flux
             ! -----------------
             f(i,j,k,MRHO) = (b1*fl1 - b2*fr1)/db12 &
                  + b1b2/db12 * drho &
                  - gg * b1b2/db12 / (1.d0 + mn) * rbdq
             f(i,j,k,MVX)  = (b1*fl2 - b2*fr2)/db12 &
                  + b1b2/db12 * drhou &
                  - gg * b1b2/db12 / (1.d0 + mn) * rbdq * ub
             f(i,j,k,MVY)  = (b1*fl3 - b2*fr3)/db12 &
                  + b1b2/db12 * drhov &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * vb + rhob * dv)
             f(i,j,k,MVZ)  = (b1*fl4 - b2*fr4)/db12 &
                  + b1b2/db12 * drhow &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * wb + rhob * dw)
             f(i,j,k,MP)  = (b1*fl5 - b2*fr5)/db12 &
                  + b1b2/db12 * drhoh &
                  - gg * b1b2/db12 / (1.d0 + mn) * (rbdq * hb + rhob * dh)
          enddo
       enddo
    enddo
  end subroutine flux
  !-----------------------------------------------------------------------
  ! find dt according to CFL condtion
  !-----------------------------------------------------------------------
  subroutine cflcond
    use grid
    integer :: i,j,k
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: vx, vy, vz, rho, p
    real(kind=DBL_KIND),dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: cs
    h = get_cellWidth()
    rho => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MRHO)
    p   => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MP)
    vx  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MVX)
    vy  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MVY)
    vz  => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MVZ)
    cs = sqrt(abs(GAMMA*p/rho))
#if defined(DIRECTIONAL_SPLIT)
    if (NDIM == 3) then
       Dtime = (CFL) / max( &
            maxval(abs(vx)+cs)/h(MX), &
            maxval(abs(vy)+cs)/h(MY), &
            maxval(abs(vz)+cs)/h(MZ) )
    elseif (NDIM == 2) then
       Dtime = (CFL) / max(  &
            maxval(abs(vx)+cs)/h(MX), &
            maxval(abs(vy)+cs)/h(MY) )
    elseif (NDIM == 1) then
       Dtime = (CFL)*minval( h(MX)/(abs(vx)+cs) )
    else
       print *, '*** error'
    endif
#elif defined(DIRECTIONAL_UNSPLIT)
    if (NDIM == 3) then
       Dtime = (CFL) / maxval( &
            (abs(vx)+cs)/h(MX) + &
            (abs(vy)+cs)/h(MY) + &
            (abs(vz)+cs)/h(MZ) )
    elseif (NDIM == 2) then
       Dtime = (CFL) / maxval(  &
            (abs(vx)+cs)/h(MX) + &
            (abs(vy)+cs)/h(MY) )
    elseif (NDIM == 1) then
       Dtime = (CFL)*minval( h(MX)/(abs(vx)+cs) )
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
    real(kind=DBL_KIND) :: dv
    dv = get_dv()

    u(:,:,:,MRHO)=q(:,:,:,MRHO)*dv
    u(:,:,:,MVX)=q(:,:,:,MVX)*u(:,:,:,MRHO)
    u(:,:,:,MVY)=q(:,:,:,MVY)*u(:,:,:,MRHO)
    u(:,:,:,MVZ)=q(:,:,:,MVZ)*u(:,:,:,MRHO)
    u(:,:,:,MP) =dv*(q(:,:,:,MP)/(GAMMA-1.d0) &
         +q(:,:,:,MRHO)*(q(:,:,:,MVX)**2+q(:,:,:,MVY)**2+q(:,:,:,MVZ)**2)/2)
  end subroutine v2u
  !-----------------------------------------------------------------------
  ! convert conservative variables u to primitive variable u
  !-----------------------------------------------------------------------
  subroutine u2v(u,q)
    use grid, only : get_dv
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: u !(IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: q !(OUT)
    real(kind=DBL_KIND) :: dv
    dv = get_dv()

    q(:,:,:,MRHO) =u(:,:,:,MRHO)/dv
    q(:,:,:,MVX) =u(:,:,:,MVX)/u(:,:,:,MRHO)
    q(:,:,:,MVY) =u(:,:,:,MVY)/u(:,:,:,MRHO)
    q(:,:,:,MVZ) =u(:,:,:,MVZ)/u(:,:,:,MRHO)
    q(:,:,:,MP) = ( u(:,:,:,MP)/dv -0.5d0*q(:,:,:,MRHO) &
         *(q(:,:,:,MVX)**2+q(:,:,:,MVY)**2+q(:,:,:,MVZ)**2) )*(GAMMA-1.0d0)
  end subroutine u2v
  !-----------------------------------------------------------------------
  ! return sound speed (c) when density (rho) is given
  !-----------------------------------------------------------------------
  subroutine get_cs(c, rho, p)
    real(kind=DBL_KIND),intent(OUT) :: c
    real(kind=DBL_KIND),intent(IN) :: rho, p
    real(kind=DBL_KIND) :: flag
    c = sqrt(abs(GAMMA*p/rho))
  end subroutine get_cs


end module flux_eos
