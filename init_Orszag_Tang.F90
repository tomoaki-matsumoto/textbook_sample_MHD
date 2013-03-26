#include "config.h"
module init
  use parameter
  implicit none
  private
  public :: init_cond
contains
  subroutine init_cond
    use grid
    use flux_eos, only : GAMMA
    integer :: i, j, k
    real(kind=DBL_KIND) :: BOXSIZE_X = PI*2.d0, BOXSIZE_Y = PI*2.d0, BOXSIZE_Z = 1.d0
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
#define RHO_(i,j,k) V(i,j,k,MRHO)
#define P_(i,j,k)   V(i,j,k,MP)
#define VX_(i,j,k)  V(i,j,k,MVX)
#define VY_(i,j,k)  V(i,j,k,MVY)
#define VZ_(i,j,k)  V(i,j,k,MVZ)
#define BX_(i,j,k)  V(i,j,k,MBX)
#define BY_(i,j,k)  V(i,j,k,MBY)
#define BZ_(i,j,k)  V(i,j,k,MBZ)
    do k = KMINGH, KMAXGH
       do j = JMINGH, JMAXGH
          do i = IMINGH, IMAXGH
             RHO_(i,j,k) = GAMMA**2
             P_(i,j,k) = GAMMA
             VX_(i,j,k) = -sin(y(j))
             VY_(i,j,k) = sin(x(i))
             VZ_(i,j,k) = 0.d0
             BX_(i,j,k) =-sin(y(j))* SQRTPI4
             BY_(i,j,k) = sin(2.d0*x(i))* SQRTPI4
             BZ_(i,j,k) = 0.d0
          end do
       end do
    end do
  end subroutine init_cond
end module init




