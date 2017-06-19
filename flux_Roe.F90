#include "config.h"
#ifndef FLUX_SCHEME_HD
ERROR: this routine is only for FLUX_SCHEME_HD. Check config.h and Makefile.
#endif
! #define WO_ENTROPY_FIX
module flux_eos
  use parameter
  implicit none
  private
!!$  real(kind=DBL_KIND),parameter :: GAMMA = 5.d0/3.d0
  real(kind=DBL_KIND),parameter :: GAMMA = 1.4d0
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
  ! get numerical flux in one dimension
  !-----------------------------------------------------------------------
  ! macro for entropy condition
#ifdef WO_ENTROPY_FIX
  ! w/o entorpy fix
#define ELMOD(EL,ELBAR,EL1, EL2) \
  EL= abs(ELBAR)
#else !WO_ENTROPY_FIX
#define ELMOD(EL,ELBAR,EL1, EL2) \
  eps=max(0.d0,((ELBAR)-(EL1)),((EL2)-(ELBAR))) ;\
  x1=0.5+sign(0.5d0,abs(ELBAR)-eps) ;\
  x2=1.d0-x1 ;\
  EL= x1*abs(ELBAR) + x2*0.5d0*((ELBAR)**2/(eps+x1)+eps)
#endif !WO_ENTROPY_FIX

  subroutine flux( ql, qr, f, ndir)
    use util
    use parameter
    real(kind=DBL_KIND),dimension(IMINGH:,JMINGH:,KMINGH:,MMIN:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(IMINGH:,JMINGH:,KMINGH:,MMIN:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(IMINGH:,JMINGH:,KMINGH:,MMIN:) :: f   ! (OUT)
    integer,intent(IN) :: ndir
    integer :: i,j,k, io, jo, ko
    real(kind=DBL_KIND) :: &
         rhol,ul,vl,wl,pl,hl,el,cl, &
         rhor,ur,vr,wr,pr,hr,er,cr, &
         dq1,dq2,dq3,dq4,dq5, &
         sql,sqr,sqa,rhob,ub,vb,wb,hb,qb2,cb2,cb,cb24i,ub2,vb2,wb2, &
         gm1,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10, &
         w1,w2,w3,w4,w5, &
         el1,el2,el3,el4,el5, &
         fl1,fl2,fl3,fl4,fl5, &
         fr1,fr2,fr3,fr4,fr5 
    real(kind=DBL_KIND) :: eps, x1, x2
    call util_arroffset(ndir,io,jo,ko)
    do k = KMIN-ko, KMAX
!$omp parallel do private(i,gm1,rhol,ul,vl,wl,pl,el,hl,cl,rhor,ur,vr,wr,pr,er,hr,cr,dq1,dq2,dq3,dq4,dq5,sql,sqr,sqa,rhob,ub,vb,wb,hb,ub2,vb2,wb2,qb2,cb2,cb,cb24i,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,w1,w2,w3,w4,w5,fl1,fl2,fl3,fl4,fl5,fr1,fr2,fr3,fr4,fr5)
       do j = JMIN-jo, JMAX
          do i = IMIN-io, IMAX
             gm1 = GAMMA - 1
             ! --------------
             ! left variables
             ! --------------
             rhol = ql(i,j,k,MRHO)
             ul = ql(i,j,k,MVX)
             vl = ql(i,j,k,MVY)
             wl = ql(i,j,k,MVZ)
             pl = ql(i,j,k,MP)
             el = rhol * (ul**2 + vl**2 + wl**2) / 2 + pl / gm1
             hl = ( el + pl )/rhol
             cl = sqrt( GAMMA * pl / rhol )
             ! --------------
             ! right variables
             ! --------------
             rhor = qr(i,j,k,MRHO)
             ur = qr(i,j,k,MVX)
             vr = qr(i,j,k,MVY)
             wr = qr(i,j,k,MVZ)
             pr = qr(i,j,k,MP)
             er = rhor * (ur**2 + vr**2 + wr**2) / 2 + pr / gm1
             hr = ( er + pr )/rhor
             cr = sqrt( GAMMA * pr / rhor )
             ! -------
             ! delta Q
             ! -------
             dq1 = rhor - rhol
             dq2 = rhor * ur - rhol * ul
             dq3 = rhor * vr - rhol * vl
             dq4 = rhor * wr - rhol * wl
             dq5 = er - el
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
             cb24i = 1/(cb2 * 4)
             ! ----------------------
             ! w = R |Lambda| R^(-1)
             ! ----------------------
             ELMOD(el1, ub-cb, ul-cl, ur-cr)
             ELMOD(el2, ub,    ul,    ur)
             ELMOD(el3, ub,    ul,    ur)
             ELMOD(el4, ub,    ul,    ur)
             ELMOD(el5, ub+cb, ul+cl, ur+cr)

             a1 = vb*dq3
             a2 = wb*dq4
             a3 = dq1*qb2
             a4 = (2*dq5 + a3 - 2*(dq2*ub + a1 + a2))
             a5 = cb*(el1 - el5)
             a6 = (el1 + el5)
             a7 = (a6 - 2*el4)
             a8 = gm1*a7*a4
             a9 = a8 - a5*(dq2 - dq1*ub)*2
             a10 = el4*qb2

             w1 = (4*cb2*dq1*el4 + a9)*cb24i

             w2 = (2*cb2*(dq2*a6 - dq1*a7*ub) & 
                  - a5*( gm1*(a3 + 2*(dq5 - a1 - a2)) &
                  + 2*( (2 - GAMMA)*dq2*ub - dq1*ub2)) &
                  + a8*ub)*cb24i

             w3 = (a9*vb + 4*cb2*(dq3*el3 + dq1*(-el3 + el4)*vb))*cb24i

             w4 = (a9*wb + 4*cb2*(dq4*el2 + dq1*(-el2 + el4)*wb))*cb24i

             w5 = (gm1*(hb*a6 - a10)*a4 &
                  + 2*cb2*(dq2*a6*ub + 2*(el3*a1 + el2*a2) &
                  + dq1*(a10 - a6*ub2 - 2*(el3*vb2 + el2*wb2))) &
                  - a5*(2*dq2*(hb - gm1*ub2) &
                  + ub*(2*dq5*gm1 - dq1*(2*hb - qb2*gm1) &
                  - 2*gm1*(a1 + a2))))*cb24i

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
             f(i,j,k,MRHO) = (fl1 + fr1 - w1)/2
             f(i,j,k,MVX)  = (fl2 + fr2 - w2)/2
             f(i,j,k,MVY)  = (fl3 + fr3 - w3)/2
             f(i,j,k,MVZ)  = (fl4 + fr4 - w4)/2
             f(i,j,k,MP)   = (fl5 + fr5 - w5)/2

          enddo
       enddo
!$omp end parallel do
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
    cs = sqrt(GAMMA*p/rho)
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
    c = sqrt(GAMMA*p/rho)
  end subroutine get_cs


end module flux_eos
