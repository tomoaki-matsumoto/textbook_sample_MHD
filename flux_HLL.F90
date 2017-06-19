#include "config.h"
#ifndef FLUX_SCHEME_HD
ERROR: this routine is only for FLUX_SCHEME_HD. Check config.h and Makefile.
#endif
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
  ! SWD(A, B) = 1.d0  if A >= B. otherwise 0.d0
#define SWD(A, B)  (0.5d0 + 0.5d0*sign(1.d0, (A)-(B)))
  ! SWR(A, B) = 1.d0 if A * B < 0, otherwise 0.d0
#define SWR(A, B)  (0.5d0 - 0.5d0*sign(1.d0, (A))*sign(1.d0, (B)))

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
         el1,el2,el3,el4,el5, &
         fl1,fl2,fl3,fl4,fl5, &
         fr1,fr2,fr3,fr4,fr5 ,&
         fa1,fa2,fa3,fa4,fa5 ,&
         ul1,ul2,ul3,ul4,ul5, &
         ur1,ur2,ur3,ur4,ur5 ,&
         ua1,ua2,ua3,ua4,ua5 ,&
         csl, csr, gm1, sr, sl
    call util_arroffset(ndir,io,jo,ko)
    do k = KMIN-ko, KMAX
!$omp parallel do private(i,gm1,rhol,ul,vl,wl,pl,el,hl,cl,rhor,ur,vr,wr,pr,er,hr,cr,ul1,ul2,ul3,ul4,ul5,fl1,fl2,fl3,fl4,fl5,ur1,ur2,ur3,ur4,ur5,fr1,fr2,fr3,fr4,fr5,sl,sr)
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
             ! UL
             ul1 = rhol
             ul2 = rhol*ul
             ul3 = rhol*vl
             ul4 = rhol*wl
             ul5 = el
             ! FL
             fl1 = rhol*ul
             fl2 = fl1*ul + pl
             fl3 = fl1*vl
             fl4 = fl1*wl
             fl5 = fl1*hl
             ! UR
             ur1 = rhor
             ur2 = rhor*ur
             ur3 = rhor*vr
             ur4 = rhor*wr
             ur5 = er
             ! FR
             fr1 = rhor*ur
             fr2 = fr1*ur + pr
             fr3 = fr1*vr
             fr4 = fr1*wr
             fr5 = fr1*hr
             ! speed
             sl = min(ul-cl, ur-cr, 0.d0)
             sr = max(ul+cl, ur+cr, 0.d0)
             ! -----------------
             ! flux
             ! -----------------
             f(i,j,k,MRHO) = (sr*fl1 - sl*fr1 + sr*sl*(ur1 - ul1))/(sr-sl)
             f(i,j,k,MVX)  = (sr*fl2 - sl*fr2 + sr*sl*(ur2 - ul2))/(sr-sl)
             f(i,j,k,MVY)  = (sr*fl3 - sl*fr3 + sr*sl*(ur3 - ul3))/(sr-sl)
             f(i,j,k,MVZ)  = (sr*fl4 - sl*fr4 + sr*sl*(ur4 - ul4))/(sr-sl)
             f(i,j,k,MP)   = (sr*fl5 - sl*fr5 + sr*sl*(ur5 - ul5))/(sr-sl)

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
