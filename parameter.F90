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
  real(kind=DBL_KIND),parameter :: CFL = 0.9d0

  integer,parameter :: NDIM = 2 ! number of dimension
  integer,parameter :: NX = 128, NY = 128, NZ = 1
  integer,parameter :: IMAX = NX-1, JMAX = NY-1, KMAX = NZ-1
  integer,parameter :: IMIN = 0,   JMIN = 0, KMIN = 0, MMIN = 0
  integer,parameter :: NGH = 2  ! ghostcell
  integer,parameter :: IMINGH = IMIN-NGH, JMINGH = JMIN-NGH, KMINGH = KMIN-NGH
  integer,parameter :: IMAXGH = IMAX+NGH, JMAXGH = JMAX+NGH, KMAXGH = KMAX+NGH
  integer,parameter :: MX = 0, MY = 1, MZ = 2
  ! for scalar advection
#if _FLUX_SCHEME_ == _SCALAR_ADVECTION_
  integer,parameter :: MRHO=0, MMAX = 0
#endif !_SCALAR_ADVECTION_
  ! for HD
#if _FLUX_SCHEME_ == _HD_
  integer,parameter :: MRHO=0, MVX=1, MVY=2, MVZ=3, MP=4, MMAX = 4
#endif !_HD_
  ! for MHD
#if _FLUX_SCHEME_ == _MHD_
  integer,parameter :: MRHO=0, MVX=1, MVY=2, MVZ=3, MBX=4, MBY=5, MBZ=6, MP=7, MDB=8, MMAX = 8
#endif !_MHD_
end module parameter
