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
  subroutine boundary_fix ( u )
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: u
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),parameter :: RHO0=1.4d0, P0=1.d0, VX0=3.d0
    integer :: istepmin, istepmax, jstepmin,jstepmax, i, j, k

    call cure_crash_v( u )

    ! inflow boundary at x=0
    u(IMINGH:IMINGH+1,:,:,MRHO) = RHO0
    u(IMINGH:IMINGH+1,:,:,MP)  = P0
    u(IMINGH:IMINGH+1,:,:,MVX) = VX0
    u(IMINGH:IMINGH+1,:,:,MVY) = 0.d0
    u(IMINGH:IMINGH+1,:,:,MVZ) = 0.d0
    ! free boundary at x=3
    do i = IMAX+1,IMAXGH
       u(i,:,:,:) = u(IMAX,:,:,:)
    enddo
    ! free boundary at jmin
    do j = JMINGH,JMIN-1
       u(:,j,:,:) = u(:,JMIN,:,:)
    enddo
    ! free boundary at jmax
    do j = JMAX+1,JMAXGH
       u(:,j,:,:) = u(:,JMAX,:,:)
    enddo
    ! z-direction is periodic (dummy)
    u(:,:,KMINGH:KMINGH+1,:) = u(:,:,KMAX-1:KMAX,:)
    u(:,:,KMAXGH-1:KMAXGH,:) = u(:,:,KMIN:KMIN+1,:)

    ! for step x in [0.62], y in [0.4,0.6]
    h = get_cellWidth()
    istepmin = 0.6d0/h(MX)
    istepmax = 0.8d0/h(MX)
    jstepmin = 0.4d0/h(MY)
    jstepmax = 0.6d0/h(MY)

    ! reflection boundary at y=0.2
#define ISTEPS istepmin:istepmax
#define JSTEPS jstepmin:jstepmax
    u(ISTEPS,JSTEPS,:,MRHO) = RHO0
    u(ISTEPS,JSTEPS,:,MP) =   P0
    u(ISTEPS,JSTEPS,:,MVX) = 0.d0
    u(ISTEPS,JSTEPS,:,MVY) = 0.d0
    u(ISTEPS,JSTEPS,:,MVZ) = 0.d0

    ! left imin
    do i = istepmin, istepmin+1
       u(i,JSTEPS,:,MRHO) = u(2*istepmin-1-i,JSTEPS,:,MRHO)
       u(i,JSTEPS,:,MP)  =  u(2*istepmin-1-i,JSTEPS,:,MP)
       u(i,JSTEPS,:,MVX) = -u(2*istepmin-1-i,JSTEPS,:,MVX)
       u(i,JSTEPS,:,MVY) =  u(2*istepmin-1-i,JSTEPS,:,MVY)
       u(i,JSTEPS,:,MVZ) =  u(2*istepmin-1-i,JSTEPS,:,MVZ)
    enddo
    ! right imax
    do i = istepmax-1, istepmax
       u(i,JSTEPS,:,MRHO) = u(2*istepmax+1-i,JSTEPS,:,MRHO)
       u(i,JSTEPS,:,MP)  =  u(2*istepmax+1-i,JSTEPS,:,MP)
       u(i,JSTEPS,:,MVX) = -u(2*istepmax+1-i,JSTEPS,:,MVX)
       u(i,JSTEPS,:,MVY) =  u(2*istepmax+1-i,JSTEPS,:,MVY)
       u(i,JSTEPS,:,MVZ) =  u(2*istepmax+1-i,JSTEPS,:,MVZ)
    enddo
    ! bottom jmin
    do j = jstepmin, jstepmin+1
       u(ISTEPS,j,:,MRHO) = u(ISTEPS,2*jstepmin-1-j,:,MRHO)
       u(ISTEPS,j,:,MP)  =  u(ISTEPS,2*jstepmin-1-j,:,MP)
       u(ISTEPS,j,:,MVX) =  u(ISTEPS,2*jstepmin-1-j,:,MVX)
       u(ISTEPS,j,:,MVY) = -u(ISTEPS,2*jstepmin-1-j,:,MVY)
       u(ISTEPS,j,:,MVZ) =  u(ISTEPS,2*jstepmin-1-j,:,MVZ)
    enddo
    ! top jmax
    do j = jstepmax-1, jstepmax
       u(ISTEPS,j,:,MRHO) = u(ISTEPS,2*jstepmax+1-j,:,MRHO)
       u(ISTEPS,j,:,MP)  =  u(ISTEPS,2*jstepmax+1-j,:,MP)
       u(ISTEPS,j,:,MVX) =  u(ISTEPS,2*jstepmax+1-j,:,MVX)
       u(ISTEPS,j,:,MVY) = -u(ISTEPS,2*jstepmax+1-j,:,MVY)
       u(ISTEPS,j,:,MVZ) =  u(ISTEPS,2*jstepmax+1-j,:,MVZ)
    enddo
  end subroutine boundary_fix
end module boundary
