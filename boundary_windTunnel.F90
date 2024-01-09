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
  subroutine boundary_fix ( u, ndir )
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: u
    integer,optional :: ndir
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),parameter :: RHO0=1.4d0, P0=1.d0, VX0=3.d0
    integer :: istepmin, istepmax, jstepmin,jstepmax, i, j, k
    logical :: bool_mx, bool_my, bool_mz
    call cure_crash_v( u )

    ! inflow boundary at x=0
    u(IMINGH:IMINGH+1,:,:,MRHO) = RHO0
    u(IMINGH:IMINGH+1,:,:,MP)  = P0
    u(IMINGH:IMINGH+1,:,:,MVX) = VX0
    u(IMINGH:IMINGH+1,:,:,MVY) = 0.d0
    u(IMINGH:IMINGH+1,:,:,MVZ) = 0.d0
    ! outgoing bondary at x=3
    do i = IMAX+1,IMAXGH
       u(i,:,:,:) = u(IMAX,:,:,:)
    enddo

    ! reflection boundary at y=0
    do j = JMINGH, JMINGH+1
       u(:,j,:,MRHO) = u(:,2*JMIN-1-j,:,MRHO)
       u(:,j,:,MP)  =  u(:,2*JMIN-1-j,:,MP)
       u(:,j,:,MVX) =  u(:,2*JMIN-1-j,:,MVX)
       u(:,j,:,MVY) = -u(:,2*JMIN-1-j,:,MVY)
       u(:,j,:,MVZ) =  u(:,2*JMIN-1-j,:,MVZ)
    enddo
    ! reflection boundary at y=1
    do j = JMAXGH-1, JMAXGH
       u(:,j,:,MRHO) = u(:,2*JMAX+1-j,:,MRHO)
       u(:,j,:,MP)  =  u(:,2*JMAX+1-j,:,MP)
       u(:,j,:,MVX) =  u(:,2*JMAX+1-j,:,MVX)
       u(:,j,:,MVY) = -u(:,2*JMAX+1-j,:,MVY)
       u(:,j,:,MVZ) =  u(:,2*JMAX+1-j,:,MVZ)
    enddo

    ! for step x>0.6, y < 0.2
    h = get_cellWidth()
    istepmin = 0.6d0/h(MX)
    istepmax = IMAXGH
    jstepmin = JMINGH
    jstepmax = 0.2d0/h(MY)

#define ISTEPS istepmin:istepmax
#define JSTEPS jstepmin:jstepmax
    u(ISTEPS,JSTEPS,:,MRHO) = RHO0
    u(ISTEPS,JSTEPS,:,MP) =   P0
    u(ISTEPS,JSTEPS,:,MVX) = 0.d0
    u(ISTEPS,JSTEPS,:,MVY) = 0.d0
    u(ISTEPS,JSTEPS,:,MVZ) = 0.d0

    ! reflection boundary at x=0.6
    do i = istepmin, istepmin+1
       u(i,JSTEPS,:,MRHO) = u(2*istepmin-1-i,JSTEPS,:,MRHO)
       u(i,JSTEPS,:,MP)  =  u(2*istepmin-1-i,JSTEPS,:,MP)
       u(i,JSTEPS,:,MVX) = -u(2*istepmin-1-i,JSTEPS,:,MVX)
       u(i,JSTEPS,:,MVY) =  u(2*istepmin-1-i,JSTEPS,:,MVY)
       u(i,JSTEPS,:,MVZ) =  u(2*istepmin-1-i,JSTEPS,:,MVZ)
    enddo

    ! reflection boundary at y=0.2
    do j = jstepmax-1, jstepmax
       u(ISTEPS,j,:,MRHO) = u(ISTEPS,2*jstepmax+1-j,:,MRHO)
       u(ISTEPS,j,:,MP)  =  u(ISTEPS,2*jstepmax+1-j,:,MP)
       u(ISTEPS,j,:,MVX) =  u(ISTEPS,2*jstepmax+1-j,:,MVX)
       u(ISTEPS,j,:,MVY) = -u(ISTEPS,2*jstepmax+1-j,:,MVY)
       u(ISTEPS,j,:,MVZ) =  u(ISTEPS,2*jstepmax+1-j,:,MVZ)
    enddo

    !corner
    u(istepmin,jstepmax,:,MRHO) = ( u(istepmin-1,jstepmax,:,MRHO) + u(istepmin,jstepmax+1,:,MRHO))/2.d0
    u(istepmin,jstepmax,:,MP  ) = ( u(istepmin-1,jstepmax,:,MP  ) + u(istepmin,jstepmax+1,:,MP  ))/2.d0
    u(istepmin,jstepmax,:,MVX ) = (-u(istepmin-1,jstepmax,:,MVX ) + u(istepmin,jstepmax+1,:,MVX ))/2.d0
    u(istepmin,jstepmax,:,MVY ) = ( u(istepmin-1,jstepmax,:,MVY ) - u(istepmin,jstepmax+1,:,MVY ))/2.d0

    u(istepmin,jstepmax-1,:,MRHO) =  u(istepmin-1,jstepmax-1,:,MRHO)
    u(istepmin,jstepmax-1,:,MP  ) =  u(istepmin-1,jstepmax-1,:,MP  )
    u(istepmin,jstepmax-1,:,MVX ) = -u(istepmin-1,jstepmax-1,:,MVX )
    u(istepmin,jstepmax-1,:,MVY ) =  u(istepmin-1,jstepmax-1,:,MVY )

    u(istepmin+1,jstepmax,:,MRHO) =  u(istepmin+1,jstepmax+1,:,MRHO)
    u(istepmin+1,jstepmax,:,MP  ) =  u(istepmin+1,jstepmax+1,:,MP  )
    u(istepmin+1,jstepmax,:,MVX ) =  u(istepmin+1,jstepmax+1,:,MVX )
    u(istepmin+1,jstepmax,:,MVY ) = -u(istepmin+1,jstepmax+1,:,MVY )

    u(istepmin+1,jstepmax-1,:,MRHO) = ( u(istepmin-2,jstepmax-1,:,MRHO) + u(istepmin+1,jstepmax+2,:,MRHO))/2.d0
    u(istepmin+1,jstepmax-1,:,MP  ) = ( u(istepmin-2,jstepmax-1,:,MP  ) + u(istepmin+1,jstepmax+2,:,MP  ))/2.d0
    u(istepmin+1,jstepmax-1,:,MVX ) = (-u(istepmin-2,jstepmax-1,:,MVX ) + u(istepmin+1,jstepmax+2,:,MVX ))/2.d0
    u(istepmin+1,jstepmax-1,:,MVY ) = ( u(istepmin-2,jstepmax-1,:,MVY ) - u(istepmin+1,jstepmax+2,:,MVY ))/2.d0


    ! z-direction is periodic (dummy)
    u(:,:,KMINGH:KMINGH+1,:) = u(:,:,KMAX-1:KMAX,:)
    u(:,:,KMAXGH-1:KMAXGH,:) = u(:,:,KMIN:KMIN+1,:)

  end subroutine boundary_fix
end module boundary
