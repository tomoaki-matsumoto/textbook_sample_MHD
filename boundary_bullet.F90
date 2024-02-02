! ---------------------------------------------------------------------------
! boundary condition
! ---------------------------------------------------------------------------
#include "config.h"
#ifndef FLUX_SCHEME_HD
ERROR: this initial condition is only for FLUX_SCHEME_HD. Check config.h and Makefile.
#endif
module boundary
  use parameter
  use grid
  use cure_crash
  implicit none
  private
  public :: boundary_fix
contains
  subroutine boundary_fix ( v )
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: v
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),parameter :: RHO0=1.4d0, P0=1.d0, VX0=3.d0
    integer :: istepmin, istepmax, jstepmin,jstepmax, i, j, k

    call cure_crash_v( v )

    ! inflow boundary at x=0
    v(IMINGH:IMINGH+1,:,:,MRHO) = RHO0
    v(IMINGH:IMINGH+1,:,:,MP)  = P0
    v(IMINGH:IMINGH+1,:,:,MVX) = VX0
    v(IMINGH:IMINGH+1,:,:,MVY) = 0.d0
    v(IMINGH:IMINGH+1,:,:,MVZ) = 0.d0
    ! free boundary at x=3
    do i = IMAX+1,IMAXGH
       v(i,:,:,:) = v(IMAX,:,:,:)
    enddo
    ! free boundary at jmin
    do j = JMINGH,JMIN-1
       v(:,j,:,:) = v(:,JMIN,:,:)
    enddo
    ! free boundary at jmax
    do j = JMAX+1,JMAXGH
       v(:,j,:,:) = v(:,JMAX,:,:)
    enddo
    ! z-direction is periodic (dummy)
    v(:,:,KMINGH:KMINGH+1,:) = v(:,:,KMAX-1:KMAX,:)
    v(:,:,KMAXGH-1:KMAXGH,:) = v(:,:,KMIN:KMIN+1,:)

    ! for step x in [0.62], y in [0.4,0.6]
    h = get_cellWidth()
    istepmin = 0.6d0/h(MX)
    istepmax = 0.8d0/h(MX)
    jstepmin = 0.4d0/h(MY)
    jstepmax = 0.6d0/h(MY)

    ! reflection boundary at y=0.2
#define ISTEPS istepmin:istepmax
#define JSTEPS jstepmin:jstepmax
    v(ISTEPS,JSTEPS,:,MRHO) = RHO0
    v(ISTEPS,JSTEPS,:,MP) =   P0
    v(ISTEPS,JSTEPS,:,MVX) = 0.d0
    v(ISTEPS,JSTEPS,:,MVY) = 0.d0
    v(ISTEPS,JSTEPS,:,MVZ) = 0.d0

    ! left imin
    do i = istepmin, istepmin+1
       v(i,JSTEPS,:,MRHO) = v(2*istepmin-1-i,JSTEPS,:,MRHO)
       v(i,JSTEPS,:,MP)  =  v(2*istepmin-1-i,JSTEPS,:,MP)
       v(i,JSTEPS,:,MVX) = -v(2*istepmin-1-i,JSTEPS,:,MVX)
       v(i,JSTEPS,:,MVY) =  v(2*istepmin-1-i,JSTEPS,:,MVY)
       v(i,JSTEPS,:,MVZ) =  v(2*istepmin-1-i,JSTEPS,:,MVZ)
    enddo
    ! right imax
    do i = istepmax-1, istepmax
       v(i,JSTEPS,:,MRHO) = v(2*istepmax+1-i,JSTEPS,:,MRHO)
       v(i,JSTEPS,:,MP)  =  v(2*istepmax+1-i,JSTEPS,:,MP)
       v(i,JSTEPS,:,MVX) = -v(2*istepmax+1-i,JSTEPS,:,MVX)
       v(i,JSTEPS,:,MVY) =  v(2*istepmax+1-i,JSTEPS,:,MVY)
       v(i,JSTEPS,:,MVZ) =  v(2*istepmax+1-i,JSTEPS,:,MVZ)
    enddo
    ! bottom jmin
    do j = jstepmin, jstepmin+1
       v(ISTEPS,j,:,MRHO) = v(ISTEPS,2*jstepmin-1-j,:,MRHO)
       v(ISTEPS,j,:,MP)  =  v(ISTEPS,2*jstepmin-1-j,:,MP)
       v(ISTEPS,j,:,MVX) =  v(ISTEPS,2*jstepmin-1-j,:,MVX)
       v(ISTEPS,j,:,MVY) = -v(ISTEPS,2*jstepmin-1-j,:,MVY)
       v(ISTEPS,j,:,MVZ) =  v(ISTEPS,2*jstepmin-1-j,:,MVZ)
    enddo
    ! top jmax
    do j = jstepmax-1, jstepmax
       v(ISTEPS,j,:,MRHO) = v(ISTEPS,2*jstepmax+1-j,:,MRHO)
       v(ISTEPS,j,:,MP)  =  v(ISTEPS,2*jstepmax+1-j,:,MP)
       v(ISTEPS,j,:,MVX) =  v(ISTEPS,2*jstepmax+1-j,:,MVX)
       v(ISTEPS,j,:,MVY) = -v(ISTEPS,2*jstepmax+1-j,:,MVY)
       v(ISTEPS,j,:,MVZ) =  v(ISTEPS,2*jstepmax+1-j,:,MVZ)
    enddo
  end subroutine boundary_fix
end module boundary
