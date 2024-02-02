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
    integer :: i, j, k
    do i = IMINGH,IMIN-1
       v(i,:,:,:) = v(IMIN,:,:,:)
    enddo
    do i = IMAX+1,IMAXGH
       v(i,:,:,:) = v(IMAX,:,:,:)
    enddo
    do j = JMINGH,JMIN-1
       v(:,j,:,:) = v(:,JMIN,:,:)
    enddo
    do j = JMAX+1,JMAXGH
       v(:,j,:,:) = v(:,JMAX,:,:)
    enddo
    do k = KMINGH,KMIN-1
       v(:,:,k,:) = v(:,:,KMIN,:)
    enddo
    do k = KMAX+1,KMAXGH
       v(:,:,k,:) = v(:,:,KMAX,:)
    enddo
  end subroutine boundary_fix
end module boundary
