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
  subroutine boundary_fix ( v )
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),intent(INOUT) :: v
    v(IMINGH:IMINGH+1,:,:,:) = v(IMAX-1:IMAX,:,:,:)
    v(IMAXGH-1:IMAXGH,:,:,:) = v(IMIN:IMIN+1,:,:,:)
    v(:,JMINGH:JMINGH+1,:,:) = v(:,JMAX-1:JMAX,:,:)
    v(:,JMAXGH-1:JMAXGH,:,:) = v(:,JMIN:JMIN+1,:,:)
    v(:,:,KMINGH:KMINGH+1,:) = v(:,:,KMAX-1:KMAX,:)
    v(:,:,KMAXGH-1:KMAXGH,:) = v(:,:,KMIN:KMIN+1,:)
  end subroutine boundary_fix
end module boundary
