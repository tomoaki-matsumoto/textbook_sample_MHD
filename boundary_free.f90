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
    integer :: i, j, k
    do i = IMINGH,IMIN-1
       u(i,:,:,:) = u(IMIN,:,:,:)
    enddo
    do i = IMAX+1,IMAXGH
       u(i,:,:,:) = u(IMAX,:,:,:)
    enddo
    do j = JMINGH,JMIN-1
       u(:,j,:,:) = u(:,JMIN,:,:)
    enddo
    do j = JMAX+1,JMAXGH
       u(:,j,:,:) = u(:,JMAX,:,:)
    enddo
    do k = KMINGH,KMIN-1
       u(:,:,k,:) = u(:,:,KMIN,:)
    enddo
    do k = KMAX+1,KMAXGH
       u(:,:,k,:) = u(:,:,KMAX,:)
    enddo
  end subroutine boundary_fix
end module boundary
