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
  subroutine boundary_fix ( v, ndir )
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: v
    integer,optional :: ndir
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    real(kind=DBL_KIND),parameter :: RHO0=1.4d0, P0=1.d0, VX0=3.d0
    integer :: istepmin, istepmax, jstepmin,jstepmax, i, j, k
    logical :: bool_mx, bool_my, bool_mz
    call cure_crash_v( v )

    ! inflow boundary at x=0
    v(IMINGH:IMINGH+1,:,:,MRHO) = RHO0
    v(IMINGH:IMINGH+1,:,:,MP)  = P0
    v(IMINGH:IMINGH+1,:,:,MVX) = VX0
    v(IMINGH:IMINGH+1,:,:,MVY) = 0.d0
    v(IMINGH:IMINGH+1,:,:,MVZ) = 0.d0
    ! outgoing bondary at x=3
    do i = IMAX+1,IMAXGH
       v(i,:,:,:) = v(IMAX,:,:,:)
    enddo

    ! reflection boundary at y=0
    do j = JMINGH, JMINGH+1
       v(:,j,:,MRHO) = v(:,2*JMIN-1-j,:,MRHO)
       v(:,j,:,MP)  =  v(:,2*JMIN-1-j,:,MP)
       v(:,j,:,MVX) =  v(:,2*JMIN-1-j,:,MVX)
       v(:,j,:,MVY) = -v(:,2*JMIN-1-j,:,MVY)
       v(:,j,:,MVZ) =  v(:,2*JMIN-1-j,:,MVZ)
    enddo
    ! reflection boundary at y=1
    do j = JMAXGH-1, JMAXGH
       v(:,j,:,MRHO) = v(:,2*JMAX+1-j,:,MRHO)
       v(:,j,:,MP)  =  v(:,2*JMAX+1-j,:,MP)
       v(:,j,:,MVX) =  v(:,2*JMAX+1-j,:,MVX)
       v(:,j,:,MVY) = -v(:,2*JMAX+1-j,:,MVY)
       v(:,j,:,MVZ) =  v(:,2*JMAX+1-j,:,MVZ)
    enddo

    ! for step x>0.6, y < 0.2
    h = get_cellWidth()
    istepmin = 0.6d0/h(MX)
    istepmax = IMAXGH
    jstepmin = JMINGH
    jstepmax = 0.2d0/h(MY)

#define ISTEPS istepmin:istepmax
#define JSTEPS jstepmin:jstepmax
    v(ISTEPS,JSTEPS,:,MRHO) = RHO0
    v(ISTEPS,JSTEPS,:,MP) =   P0
    v(ISTEPS,JSTEPS,:,MVX) = 0.d0
    v(ISTEPS,JSTEPS,:,MVY) = 0.d0
    v(ISTEPS,JSTEPS,:,MVZ) = 0.d0

    ! reflection boundary at x=0.6
    do i = istepmin, istepmin+1
       v(i,JSTEPS,:,MRHO) = v(2*istepmin-1-i,JSTEPS,:,MRHO)
       v(i,JSTEPS,:,MP)  =  v(2*istepmin-1-i,JSTEPS,:,MP)
       v(i,JSTEPS,:,MVX) = -v(2*istepmin-1-i,JSTEPS,:,MVX)
       v(i,JSTEPS,:,MVY) =  v(2*istepmin-1-i,JSTEPS,:,MVY)
       v(i,JSTEPS,:,MVZ) =  v(2*istepmin-1-i,JSTEPS,:,MVZ)
    enddo

    ! reflection boundary at y=0.2
    do j = jstepmax-1, jstepmax
       v(ISTEPS,j,:,MRHO) = v(ISTEPS,2*jstepmax+1-j,:,MRHO)
       v(ISTEPS,j,:,MP)  =  v(ISTEPS,2*jstepmax+1-j,:,MP)
       v(ISTEPS,j,:,MVX) =  v(ISTEPS,2*jstepmax+1-j,:,MVX)
       v(ISTEPS,j,:,MVY) = -v(ISTEPS,2*jstepmax+1-j,:,MVY)
       v(ISTEPS,j,:,MVZ) =  v(ISTEPS,2*jstepmax+1-j,:,MVZ)
    enddo

    !corner
    v(istepmin,jstepmax,:,MRHO) = ( v(istepmin-1,jstepmax,:,MRHO) + v(istepmin,jstepmax+1,:,MRHO))/2.d0
    v(istepmin,jstepmax,:,MP  ) = ( v(istepmin-1,jstepmax,:,MP  ) + v(istepmin,jstepmax+1,:,MP  ))/2.d0
    v(istepmin,jstepmax,:,MVX ) = (-v(istepmin-1,jstepmax,:,MVX ) + v(istepmin,jstepmax+1,:,MVX ))/2.d0
    v(istepmin,jstepmax,:,MVY ) = ( v(istepmin-1,jstepmax,:,MVY ) - v(istepmin,jstepmax+1,:,MVY ))/2.d0

    v(istepmin,jstepmax-1,:,MRHO) =  v(istepmin-1,jstepmax-1,:,MRHO)
    v(istepmin,jstepmax-1,:,MP  ) =  v(istepmin-1,jstepmax-1,:,MP  )
    v(istepmin,jstepmax-1,:,MVX ) = -v(istepmin-1,jstepmax-1,:,MVX )
    v(istepmin,jstepmax-1,:,MVY ) =  v(istepmin-1,jstepmax-1,:,MVY )

    v(istepmin+1,jstepmax,:,MRHO) =  v(istepmin+1,jstepmax+1,:,MRHO)
    v(istepmin+1,jstepmax,:,MP  ) =  v(istepmin+1,jstepmax+1,:,MP  )
    v(istepmin+1,jstepmax,:,MVX ) =  v(istepmin+1,jstepmax+1,:,MVX )
    v(istepmin+1,jstepmax,:,MVY ) = -v(istepmin+1,jstepmax+1,:,MVY )

    v(istepmin+1,jstepmax-1,:,MRHO) = ( v(istepmin-2,jstepmax-1,:,MRHO) + v(istepmin+1,jstepmax+2,:,MRHO))/2.d0
    v(istepmin+1,jstepmax-1,:,MP  ) = ( v(istepmin-2,jstepmax-1,:,MP  ) + v(istepmin+1,jstepmax+2,:,MP  ))/2.d0
    v(istepmin+1,jstepmax-1,:,MVX ) = (-v(istepmin-2,jstepmax-1,:,MVX ) + v(istepmin+1,jstepmax+2,:,MVX ))/2.d0
    v(istepmin+1,jstepmax-1,:,MVY ) = ( v(istepmin-2,jstepmax-1,:,MVY ) - v(istepmin+1,jstepmax+2,:,MVY ))/2.d0


    ! z-direction is periodic (dummy)
    v(:,:,KMINGH:KMINGH+1,:) = v(:,:,KMAX-1:KMAX,:)
    v(:,:,KMAXGH-1:KMAXGH,:) = v(:,:,KMIN:KMIN+1,:)

  end subroutine boundary_fix
end module boundary
