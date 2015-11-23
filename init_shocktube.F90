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
    real(kind=DBL_KIND) :: BOXSIZE_X = 5.d0, BOXSIZE_Y = 1.d0, BOXSIZE_Z = 1.d0
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

#ifdef FLUX_SCHEME_SCALAR_ADVECTION
             if ( x(i) <= 0.d0 ) then
                V(i,j,k,MRHO) = 2.d0
             else
                V(i,j,k,MRHO) = 1.d0
             endif
#endif !FLUX_SCHEME_SCALAR_ADVECTION


#ifdef FLUX_SCHEME_HD
             if ( x(i) <= 0.d0 ) then
                V(i,j,k,MRHO) = 1.d0
!                V(i,j,k,MVX) = 0.9d0  ! for entropy fix
                V(i,j,k,MP) = 1.d0
             else
                V(i,j,k,MRHO) = 0.125d0
!                V(i,j,k,MVX) = 0.9d0   ! for entropy fix
                V(i,j,k,MP) = 0.1d0
             end if
#endif !FLUX_SCHEME_HD


#ifdef FLUX_SCHEME_MHD
             if ( x(i) <= 0.d0 ) then
                V(i,j,k,MRHO) = 1.d0
                V(i,j,k,MP) = 1.d0
                V(i,j,k,MBX) = 0.75d0 * SQRTPI4
                V(i,j,k,MBY) = 1.d0 * SQRTPI4
                V(i,j,k,MBZ) = 0.d0
             else
                V(i,j,k,MRHO) = 0.125d0
                V(i,j,k,MP) = 0.1d0
                V(i,j,k,MBX) = 0.75d0 * SQRTPI4
                V(i,j,k,MBY) = -1.d0 * SQRTPI4
                V(i,j,k,MBZ) = 0.d0
             end if
#endif !FLUX_SCHEME_MHD

          end do
       end do
    end do
  end subroutine init_cond
end module init




