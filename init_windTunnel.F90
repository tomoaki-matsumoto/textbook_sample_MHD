! initial condition for wind tunnel problem with a step.
#include "config.h"
#ifndef FLUX_SCHEME_HD
ERROR: this initial condition is only for FLUX_SCHEME_HD. Check config.h and Makefile.
#endif
module init
  use parameter
  implicit none
  private
  public :: init_cond
contains
  subroutine init_cond
    use grid
    integer :: i, j, k
    real(kind=DBL_KIND) :: BOXSIZE_X = 3.d0, BOXSIZE_Y = 1.d0, BOXSIZE_Z = 1.d0
    real(kind=DBL_KIND) :: dx, dy, dz
    real(kind=DBL_KIND) :: bx, by, bz
    real(kind=DBL_KIND) :: ic0, jc0, kc0
    integer :: istepmin, istepmax, jstepmin, jstepmax
    ! define parameters
    Time = 0.d0
    Step = 0
    Dtime = 0.1d0
    ! define coordinates
    dx = BOXSIZE_X/(IMAX-IMIN+1)
    dy = BOXSIZE_Y/(JMAX-JMIN+1)
    dz = BOXSIZE_Z/(KMAX-KMIN+1)
    ! origin in cell number
    ic0 = -0.5d0
    jc0 = -0.5d0
    kc0 = -0.5d0
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
    V(:,:,:,MRHO) = 1.4d0
    V(:,:,:,MP)   = 1.d0
    V(:,:,:,MVX)  = 3.d0
    V(:,:,:,MVY)  = 0.d0
    V(:,:,:,MVZ)  = 0.d0

    istepmin = 0.6d0/dx
    istepmax = IMAXGH
    jstepmin = JMINGH
    jstepmax = 0.2d0/dy
#define ISTEPS istepmin:istepmax
#define JSTEPS jstepmin:jstepmax
    V(ISTEPS,JSTEPS,:,MVX) = 0.d0

  end subroutine init_cond
end module init




