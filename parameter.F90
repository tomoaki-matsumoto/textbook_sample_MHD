#include "config.h"
module parameter
  implicit none
  public
  integer,parameter :: DBL_KIND = 8
  real(kind=DBL_KIND),parameter :: PI = atan(1.d0)*4.d0
  real(kind=DBL_KIND),parameter :: PI4 = PI * 4.d0
  real(kind=DBL_KIND),parameter :: PI4I = 1.d0/PI4
  real(kind=DBL_KIND),parameter :: SQRTPI = sqrt(PI)
  real(kind=DBL_KIND),parameter :: SQRTPI4 = sqrt(PI4)
  real(kind=DBL_KIND),parameter :: CFL = 0.5d0

  ! Condition for a halt
  integer,parameter :: STEPMAX = 1000000
!!$  real(kind=DBL_KIND),parameter :: T_LAST = 0.2d0
!!$  real(kind=DBL_KIND),parameter :: T_LAST = PI
  real(kind=DBL_KIND),parameter :: T_LAST = 1.d0

#ifndef CONVERGENCE_TEST
  integer,parameter :: NDIM = 1 ! number of dimension
  integer,parameter :: NX = 64, NY = 1, NZ = 1
#endif !CONVERGENCE_TEST
  integer,parameter :: IMAX = NX-1, JMAX = NY-1, KMAX = NZ-1
  integer,parameter :: IMIN = 0,   JMIN = 0, KMIN = 0, MMIN = 0
  integer,parameter :: NGH = 2  ! ghostcell
  integer,parameter :: IMINGH = IMIN-NGH, JMINGH = JMIN-NGH, KMINGH = KMIN-NGH
  integer,parameter :: IMAXGH = IMAX+NGH, JMAXGH = JMAX+NGH, KMAXGH = KMAX+NGH
  integer,parameter :: MX = 0, MY = 1, MZ = 2
  ! for scalar advection
#if defined(FLUX_SCHEME_SCALAR_ADVECTION)
  integer,parameter :: MRHO=0, MMAX = 0
  ! for HD
#elif defined(FLUX_SCHEME_HD)
  integer,parameter :: MRHO=0, MVX=1, MVY=2, MVZ=3, MP=4, MMAX = 4
  ! for MHD
#elif defined(FLUX_SCHEME_MHD)
  integer,parameter :: MRHO=0, MVX=1, MVY=2, MVZ=3, MBX=4, MBY=5, MBZ=6, MP=7, MDB=8, MMAX = 8
#else
  ERROR
#endif
end module parameter
