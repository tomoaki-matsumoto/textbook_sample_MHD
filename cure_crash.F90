!
! cure clashing simulation
!  
#include "config.h"

! Cure by the way so that the momentum is conserved.
! #define CURE_MOMENTUM_CONSERVATION

! Cure by the minimum operations. The influence of the cure is minimum.
#define CURE_MINIMUM

module cure_crash
  implicit none
  private
#ifndef FLUX_SCHEME_SCALAR_ADVECTION
  public :: cure_crash_v
contains
  subroutine cure_crash_v( v )
    use parameter
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: v
    real(kind=DBL_KIND),parameter :: FLOOR_RHO = 0.d0, FLOOR_P = 0.d0
    integer :: i, j, k, is, ie, js, je, ks, ke, ncell
    real(kind=DBL_KIND) :: rho_ave, p_ave, vx_ave, vy_ave, vz_ave
#define RHO(i,j,k) v(i,j,k,MRHO)
#define P(i,j,k) v(i,j,k,MP)
#define VX(i,j,k) v(i,j,k,MVX)
#define VY(i,j,k) v(i,j,k,MVY)
#define VZ(i,j,k) v(i,j,k,MVZ)

#ifdef CURE_MOMENTUM_CONSERVATION
    do k=KMIN,KMAX
!$omp parallel do private(i,is,ie,js,je,ks,ke,ncell,rho_ave,p_ave,vx_ave,vy_ave,vz_ave)
       do j=JMIN,JMAX
          do i=IMIN,IMAX
             if (RHO(i,j,k) <= FLOOR_RHO .or. P(i,j,k) <= FLOOR_P) then
                write(*,"('floor at (', 3I4, '), rho, p =', 2E13.5)") i,j,k, RHO(i,j,k), P(i,j,k)
                is = max(i-1, IMIN)
                ie = min(i+1, IMAX)
                js = max(j-1, JMIN)
                je = min(j+1, JMAX)
                ks = max(k-1, KMIN)
                ke = min(k+1, KMAX)
                ncell = (ie-is+1)*(je-js+1)*(ke-ks+1)
                rho_ave = sum(RHO(is:ie,js:je,ks:ke))/ncell
                p_ave = sum(P(is:ie,js:je,ks:ke))/ncell
                vx_ave = sum(RHO(is:ie,js:je,ks:ke)*VX(is:ie,js:je,ks:ke))/ncell / rho_ave
                vy_ave = sum(RHO(is:ie,js:je,ks:ke)*VY(is:ie,js:je,ks:ke))/ncell / rho_ave
                vz_ave = sum(RHO(is:ie,js:je,ks:ke)*VZ(is:ie,js:je,ks:ke))/ncell / rho_ave
                RHO(is:ie,js:je,ks:ke) = rho_ave
                P(is:ie,js:je,ks:ke) = p_ave
                VX(is:ie,js:je,ks:ke) = vx_ave
                VY(is:ie,js:je,ks:ke) = vy_ave
                VZ(is:ie,js:je,ks:ke) = vz_ave
             endif
          end do
       end do
    end do
#endif !CURE_MOMENTUM_CONSERVATION

#ifdef CURE_MINIMUM
    do k=KMIN,KMAX
!$omp parallel do private(i,is,ie,js,je,ks,ke,ncell,rho_ave,p_ave,vx_ave,vy_ave,vz_ave)
       do j=JMIN,JMAX
          do i=IMIN,IMAX
             if (RHO(i,j,k) <= FLOOR_RHO) then
                write(*,"('floor at (', 3I4, '), rho =', 2E13.5)") i,j,k, RHO(i,j,k)
                is = max(i-1, IMIN)
                ie = min(i+1, IMAX)
                js = max(j-1, JMIN)
                je = min(j+1, JMAX)
                ks = max(k-1, KMIN)
                ke = min(k+1, KMAX)
                ncell = (ie-is+1)*(je-js+1)*(ke-ks+1)
                rho_ave = sum(RHO(is:ie,js:je,ks:ke))/ncell
                RHO(i,j,k) = rho_ave
             endif
             if (P(i,j,k) <= FLOOR_P) then
                write(*,"('floor at (', 3I4, '), p =', 2E13.5)") i,j,k, P(i,j,k)
                is = max(i-1, IMIN)
                ie = min(i+1, IMAX)
                js = max(j-1, JMIN)
                je = min(j+1, JMAX)
                ks = max(k-1, KMIN)
                ke = min(k+1, KMAX)
                ncell = (ie-is+1)*(je-js+1)*(ke-ks+1)
                p_ave = sum(P(is:ie,js:je,ks:ke))/ncell
                P(i,j,k) = p_ave
             endif
          end do
       end do
    end do
#endif !CURE_MINIMUM

  end subroutine cure_crash_v
#endif !FLUX_SCHEME_SCALAR_ADVECTION
end module cure_crash
