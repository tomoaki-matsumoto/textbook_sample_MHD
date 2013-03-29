! initial condition for advection problem for contact discontinuity
#include "config.h"
module init
  use parameter
  implicit none
  private
  public :: init_cond
contains
  subroutine init_cond
    use grid
    integer :: i, j, k
    real(kind=DBL_KIND) :: BOXSIZE_X = 1.d0, BOXSIZE_Y = 1.d0, BOXSIZE_Z = 1.d0
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
!!$    ic0 = -0.5d0
!!$    jc0 = -0.5d0
!!$    kc0 = -0.5d0
    do i = IMINGH, IMAXGH
       X(i) = (i-ic0) * dx
    end do
    do j = JMINGH, JMAXGH
       Y(j) = (j-jc0) * dy
    end do
    do k = KMINGH, KMAXGH
       Z(k) = (k-kc0) * dz
    end do
#if _FLUX_SCHEME_ == _HD_
    ! define physical variables
    V(:,:,:,MVX)  = 1.d0
    V(:,:,:,MVY)  = 1.d0
    V(:,:,:,MVZ)  = 0.d0
    V(:,:,:,MP)   = 1.d0
#endif !_HD_
    V(:,:,:,MRHO) = 1.d0
    do k = KMINGH, KMAXGH
       do j = JMINGH, JMAXGH
          do i = IMINGH, IMAXGH
             if ( x(i) >= -0.1 .and. x(i) <= 0.1 .and. &
                  y(j) >= -0.1 .and. y(j) <= 0.1 ) then
                V(i,j,k,MRHO) = 2.d0
             end if
          end do
       end do
    end do
  end subroutine init_cond
end module init




