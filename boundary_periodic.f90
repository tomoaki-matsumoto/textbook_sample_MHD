! ---------------------------------------------------------------------------
! boundary condition
! ---------------------------------------------------------------------------
module boundary
  use parameter
  use grid
  implicit none
  private
  public :: boundary_fix
contains
  subroutine boundary_fix ( u )
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: u
    u(IMINGH:IMINGH+1,:,:,:) = u(IMAX-1:IMAX,:,:,:)
    u(IMAXGH-1:IMAXGH,:,:,:) = u(IMIN:IMIN+1,:,:,:)
    u(:,JMINGH:JMINGH+1,:,:) = u(:,JMAX-1:JMAX,:,:)
    u(:,JMAXGH-1:JMAXGH,:,:) = u(:,JMIN:JMIN+1,:,:)
    u(:,:,KMINGH:KMINGH+1,:) = u(:,:,KMAX-1:KMAX,:)
    u(:,:,KMAXGH-1:KMAXGH,:) = u(:,:,KMIN:KMIN+1,:)
  end subroutine boundary_fix
end module boundary
