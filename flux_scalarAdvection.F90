#include "config.h"
#ifndef FLUX_SCHEME_SCALAR_ADVECTION
ERROR
#endif
module flux_eos
  use parameter
  implicit none
  private
  real(kind=DBL_KIND),parameter :: CS = 1.d0
  public :: cyclecomp, flux, cflcond, v2u, u2v, CS
contains
  !-----------------------------------------------------------------------
  ! rotates components of system equation
  !-----------------------------------------------------------------------
  function cyclecomp(ncrd, invert) result( mcycle )
    integer,intent(IN) :: ncrd
    integer,intent(IN),optional :: invert
    integer,dimension(MMIN:MMAX) :: mcycle
    integer :: m
    mcycle = MMIN
  end function cyclecomp
  !-----------------------------------------------------------------------
  ! get numerical flux in one dimension
  !-----------------------------------------------------------------------
  subroutine flux( ql, qr, f)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: ql  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: qr  ! (IN)
    real(kind=DBL_KIND),dimension(:,:,:,MMIN:) :: f   ! (OUT)
    if (CS <= 0) then
       f = qr
    else
       f = ql
    endif
  end subroutine flux
  !-----------------------------------------------------------------------
  ! find dt according to CFL condtion
  !-----------------------------------------------------------------------
  subroutine cflcond
    use grid
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),dimension(:,:,:),pointer :: rho
    h = get_cellWidth()
    rho => V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MRHO)
#if defined(DIRECTIONAL_SPLIT)
    if (NDIM == 3) then
       Dtime = (CFL) / CS * min( h(MX), h(MY), h(MZ) )
    elseif (NDIM == 2) then
       Dtime = (CFL) / CS * min( h(MX), h(MY) )
    elseif (NDIM == 1) then
       Dtime = (CFL) / CS * h(MX)
    else
       print *, '*** error'
    endif
#elif defined(DIRECTIONAL_UNSPLIT)
    if (NDIM == 3) then
       Dtime = (CFL) / CS /( 1.d0/h(MX) + 1.d0/h(MY) + 1.d0/h(MZ) )
    elseif (NDIM == 2) then
       Dtime = (CFL) / CS /( 1.d0/h(MX) + 1.d0/h(MY))
    elseif (NDIM == 1) then
       Dtime = (CFL) * h(MX)/ CS
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
  end subroutine u2v

end module flux_eos
