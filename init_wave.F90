#include "config.h"
! initial condition for wave propergation
module init
  use parameter
  implicit none
  private
  public :: init_cond
contains
  subroutine init_cond
    use grid
    integer :: i, j, k
    real(kind=DBL_KIND),parameter :: BOXSIZE_X = 1.d0, BOXSIZE_Y = 1.d0, BOXSIZE_Z = 1.d0
#if defined(FLUX_SCHEME_SCALAR_ADVECTION)
    real(kind=DBL_KIND),parameter :: CS = 1.d0
    real(kind=DBL_KIND),parameter :: RHO0 = 0.d0
    real(kind=DBL_KIND),parameter :: AMPLITUDE=1.d0
#elif defined(FLUX_SCHEME_HD)
    real(kind=DBL_KIND),parameter :: CS = 1.d0
    real(kind=DBL_KIND),parameter :: GAMMA = 5.d0/3.d0
    real(kind=DBL_KIND),parameter :: RHO0 = 1.d0
    real(kind=DBL_KIND),parameter :: P0 = CS**2 * RHO0 / GAMMA
    real(kind=DBL_KIND),parameter :: AMPLITUDE=1.d-6
#endif
    real(kind=DBL_KIND) :: dx, dy, dz
    real(kind=DBL_KIND) :: bx, by, bz
    real(kind=DBL_KIND) :: ic0, jc0, kc0
    ! define parameters
    Time = 0.d0
    Step = 0
    Dtime = 0.1d0
    ! define coordinates
    dx = BOXSIZE_X/(IMAX-IMIN+1)
    dy = BOXSIZE_Y/(JMAX-JMIN+1)
    dz = BOXSIZE_Z/(KMAX-KMIN+1)
    ! origin in cell number
    ic0 = (IMAX-IMIN)/2.d0
    jc0 = (JMAX-JMIN)/2.d0
    kc0 = (KMAX-KMIN)/2.d0
    do i = IMINGH, IMAXGH
       X(i) = (i-ic0) * dx
    end do
    do j = JMINGH, JMAXGH
       Y(j) = (j-jc0) * dy
    end do
    do k = KMINGH, KMAXGH
       Z(k) = (k-kc0) * dz
    end do
    ! define physical variables
    V(:,:,:,:) = 0.d0           ! initialize
    do k = KMINGH, KMAXGH
       do j = JMINGH, JMAXGH
          do i = IMINGH, IMAXGH
#if defined(FLUX_SCHEME_SCALAR_ADVECTION)
             V(i,j,k,MRHO) = sin(x(i)/BOXSIZE_X*2.d0*PI)*AMPLITUDE+RHO0
#elif defined(FLUX_SCHEME_HD)
             V(i,j,k,MRHO) = RHO0 + RHO0*AMPLITUDE*sin(x(i)/BOXSIZE_X*2.d0*PI)
             V(i,j,k,MP)   = P0 + RHO0*CS**2*AMPLITUDE*sin(x(i)/BOXSIZE_X*2.d0*PI)
             V(i,j,k,MVX)  = CS*AMPLITUDE*sin(x(i)/BOXSIZE_X*2.d0*PI)
#endif
          end do
       end do
    end do
  end subroutine init_cond
end module init




