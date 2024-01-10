! Initial condition for Kelvin-Helmholtz instability
! Reference : McNally, C. P., Lyra, W., & Passy, J.-C. 2012, APJS, 201, 18. doi:10.1088/0067-0049/201/2/18
#include "config.h"
#ifndef FLUX_SCHEME_HD
ERROR: this initial condition is only for FLUX_SCHEME_HD. Check config.h and Makefile.
#endif
module init
  use parameter
  implicit none
  private
  public :: init_cond
contains
  subroutine init_cond
    use grid
    integer :: i, j, k
    real(kind=DBL_KIND) :: BOXSIZE_X = 1.d0, BOXSIZE_Y = 1.d0, BOXSIZE_Z = 1.d0
    real(kind=DBL_KIND) :: dx, dy, dz
    real(kind=DBL_KIND) :: ic0, jc0, kc0
    real(kind=DBL_KIND) :: rho1, rho2, rhom, u1, u2, um, lsmooth
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
    rho1 = 1.d0
    rho2 = 2.d0
    rhom = (rho1 - rho2)/2.d0
    u1 = 0.5d0
    u2 = -0.5d0
    um = (u1 - u2)/2.d0
    lsmooth = 0.025d0
    do j=JMINGH, JMAXGH
       if (1./4. > y(j)) then
          RHO_(:,j,:) = rho1 - rhom * exp((y(j)-1./4.)/lsmooth)
          VX_(:,j,:)  = u1 - um * exp((y(j)-1./4.)/lsmooth)
       elseif (1./2. > y(j) .and. y(j) >= 1./4.) then
          RHO_(:,j,:) = rho2 + rhom * exp((-y(j)+1./4.)/lsmooth)
          VX_(:,j,:) = u2 + um * exp((-y(j)+1./4.)/lsmooth)
       elseif (3./4. > y(j) .and. y(j) >= 1./2.) then
          RHO_(:,j,:) = rho2 + rhom * exp(-(3./4.-y(j))/lsmooth)
          VX_(:,j,:) = u2 + um * exp(-(3./4.-y(j))/lsmooth)
       else
          RHO_(:,j,:) = rho1 - rhom * exp(-(y(j)-3./4.)/lsmooth)
          VX_(:,j,:) = u1 - um * exp(-(y(j)-3./4.)/lsmooth)
       endif
    end do
    do i = IMINGH, IMAXGH
       VY_(i,:,:) = 0.01d0 * sin(PI4 * x(i))
    end do
    VZ_(:,:,:) = 0.d0
    P_(:,:,:) = 2.5d0
  end subroutine init_cond
end module init




