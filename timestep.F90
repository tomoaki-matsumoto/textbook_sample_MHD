! Directional split : unsplit, split
! Time marching     : Euler, Predictor-corrector (PC2), Runge-Kutta (RK2), Runge-Kutta (RK3)
! Reconstruction    : NON, MUSCL2, MUSCL3
! Limiter           : MINMOD, SuperBee
#include "config.h"
module timestep
  use parameter
  implicit none
  private
  public :: step_full
contains
  ! ---------------------------------------------------------------------------
  ! 
  ! ---------------------------------------------------------------------------
  subroutine step_full

#if (_DIRECTIONAL_SPLIT_ == _UNSPLIT_) && (_TIME_MARCHING_ == _EULER_)
    call step_unsplit_euler
#endif

#if (_DIRECTIONAL_SPLIT_ == _UNSPLIT_) && (_TIME_MARCHING_ == _PC2_)
    call step_unsplit_pc2
#endif

#if (_DIRECTIONAL_SPLIT_ == _UNSPLIT_) && (_TIME_MARCHING_ == _RK2_)
    call step_unsplit_rk2
#endif

#if (_DIRECTIONAL_SPLIT_ == _UNSPLIT_) && (_TIME_MARCHING_ == _RK3_)
    call step_unsplit_rk3
#endif

#if (_DIRECTIONAL_SPLIT_ == _SPLIT_)
    call step_split
#endif

  end subroutine step_full
  ! ---------------------------------------------------------------------------
  ! 2nd order Runge-Kutta: case of Heun's method (RK2 alpha = 1)
  ! ---------------------------------------------------------------------------
  subroutine step_unsplit_rk2
    use grid
    use flux_eos
    use boundary
    ! L^n = L(U^n)
    ! U* = U^n + dt L^n
    ! L* = L(V*)
    ! U^(n+1) = 1/2 ( U^n + U* + dt L* )
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux
    W = U
    call w_update(Dtime)        ! W := U*
    call u2v(W, V)
    call boundary_fix(V)
    call get_flux
    W = (U + W)*0.5d0
    call w_update(Dtime*0.5d0)
    call u2v(W, V)
  end subroutine step_unsplit_rk2
  ! ---------------------------------------------------------------------------
  ! 3rd order Runge-Kutta (see PLUTO users guide)
  ! ---------------------------------------------------------------------------
  subroutine step_unsplit_rk3
    use grid
    use flux_eos
    use boundary
    ! U*  = U^n + dt L^n
    ! U** = 1/4 ( 3U^n + U* + dt L* )
    ! U^(n+1) = 1/3 ( U^n + 2U** + 2 dt L** )
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux
    W = U
    call w_update(Dtime)        ! W := U*, U = U^n
    call u2v(W, V)              ! V := V*
    call boundary_fix(V)
    call get_flux
    W = ( 3.d0*U + W ) *0.25d0 
    call w_update(Dtime*0.25d0) ! W := U**
    call u2v(W, V)              ! V := V**
    call boundary_fix(V) 
    call get_flux
    W = ( U + 2.d0*W ) /3.d0
    call w_update(Dtime * 2.d0/3.d0)
    call u2v(W, V)              ! V := V**
  end subroutine step_unsplit_rk3
  ! ---------------------------------------------------------------------------
  ! 2nd order method: unsplit predictor-corrector (equivalent to RK2 alpha=1/2)
  ! ---------------------------------------------------------------------------
  subroutine step_unsplit_pc2
    use grid
    use flux_eos
    use boundary
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux
    W = U
    call w_update(Dtime*0.5d0)
    call u2v(W, V)
    call boundary_fix(V)
    call get_flux
    W = U
    call w_update(Dtime)
    call u2v(W, V)
  end subroutine step_unsplit_pc2
  ! ---------------------------------------------------------------------------
  ! 1st order method: euler method
  ! ---------------------------------------------------------------------------
  subroutine step_unsplit_euler
    use grid
    use flux_eos
    use boundary
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux
    W = U
    call w_update(Dtime)
    call u2v(W, V)
  end subroutine step_unsplit_euler
  ! ---------------------------------------------------------------------------
  ! Macro for step_split_1d
  ! ---------------------------------------------------------------------------
#if _TIME_MARCHING_ == _EULER_
#define STEP_SPLIT_1D(n) step_split_euler_1d(n)
#endif
#if _TIME_MARCHING_ == _PC2_
#define STEP_SPLIT_1D(n) step_split_pc2_1d(n)
#endif
#if _TIME_MARCHING_ == _RK2_
#define STEP_SPLIT_1D(n) step_split_rk2_1d(n)
#endif
#if _TIME_MARCHING_ == _RK3_
#define STEP_SPLIT_1D(n) step_split_rk3_1d(n)
#endif
  ! ---------------------------------------------------------------------------
  ! fractional time step
  ! ---------------------------------------------------------------------------
  subroutine step_split
#ifdef _MHD_
    use grid, only: V, W, F, Dtime
    use flux_eos, only : source_b, v2u, u2v
#endif
    if (NDIM == 1) then
       call STEP_SPLIT_1D(MX)
    elseif (NDIM == 2) then
       call step_split_2d
    elseif (NDIM == 3) then
       call step_split_3d
    else
       print *, '*** error in step_fractional'
    endif
#ifdef _MHD_
    call v2u(V, W)
    call source_b(F, W, Dtime)
    call u2v(W, V)
#endif
  end subroutine step_split
  ! ---------------------------------------------------------------------------
  ! fractional time step for 3D
  ! ---------------------------------------------------------------------------
  subroutine step_split_3d
    use grid
    integer,dimension(MX:MZ),parameter :: ncycle = (/ MX, MY, MZ /)
    integer,dimension(MX:MZ) :: ndirs
    integer :: n, ns, ne, nskip, nord_R, nord_L
    ndirs = cshift( ncycle, mod(Step/2, 3))
    nord_R = mod(Step, 2)       ! right-hand order
    nord_L = 1 - nord_R         ! left-hand order
    ns = nord_R*MX + nord_L*MZ
    ne = nord_R*MZ + nord_L*MX
    nskip = nord_R - nord_L
    do n = ns, ne, nskip
       call STEP_SPLIT_1D( ndirs(n) )
    end do
  end subroutine step_split_3d
  ! ---------------------------------------------------------------------------
  ! fractional time step for 2D
  ! ---------------------------------------------------------------------------
  subroutine step_split_2d
    integer :: n
    do n = MX, MY
       call STEP_SPLIT_1D( n )
    end do
  end subroutine step_split_2d
  ! ---------------------------------------------------------------------------
  ! Runge-Kutta 2 in one dimension
  ! ---------------------------------------------------------------------------
  subroutine step_split_rk2_1d(ndir)
    use grid
    use flux_eos
    use boundary
    integer,intent(IN) :: ndir
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = U
    call w_update_ndir(Dtime, ndir)        ! W := U*
    call u2v(W, V)
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = (U + W)*0.5d0
    call w_update_ndir(Dtime*0.5d0, ndir)
    call u2v(W, V)
  end subroutine step_split_rk2_1d
  ! ---------------------------------------------------------------------------
  ! Runge-Kutta 3 in one dimension
  ! ---------------------------------------------------------------------------
  subroutine step_split_rk3_1d(ndir)
    use grid
    use flux_eos
    use boundary
    integer,intent(IN) :: ndir
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = U
    call w_update_ndir(Dtime,ndir)        ! W := U*, U = U^n
    call u2v(W, V)              ! V := V*
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = ( 3.d0*U + W ) *0.25d0 
    call w_update_ndir(Dtime*0.25d0,ndir) ! W := U**
    call u2v(W, V)              ! V := V**
    call boundary_fix(V) 
    call get_flux_ndir(ndir)
    W = ( U + 2.d0*W ) /3.d0
    call w_update_ndir(Dtime * 2.d0/3.d0, ndir)
    call u2v(W, V)              ! V := V**
  end subroutine step_split_rk3_1d
  ! ---------------------------------------------------------------------------
  ! predictor corrector in one dimension
  ! ---------------------------------------------------------------------------
  subroutine step_split_pc2_1d(ndir)
    use grid
    use flux_eos
    use boundary
    integer,intent(IN) :: ndir
    ! predictor
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = U
    call w_update_ndir(Dtime*0.5d0, ndir)
    call u2v(W, V)
    ! corrector
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = U
    call w_update_ndir(Dtime, ndir)
    call u2v(W, V)
  end subroutine step_split_pc2_1d
  ! ---------------------------------------------------------------------------
  ! Eular method in one dimension
  ! ---------------------------------------------------------------------------
  subroutine step_split_euler_1d(ndir)
    use grid
    use flux_eos
    use boundary
    integer,intent(IN) :: ndir
    ! predictor
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux_ndir(ndir)
    W = U
    call w_update_ndir(Dtime, ndir)
    call u2v(W, V)
  end subroutine step_split_euler_1d
  ! ---------------------------------------------------------------------------
  ! update u by flux for all directions
  ! ---------------------------------------------------------------------------
  subroutine w_update( dt )
    use grid
#ifdef _MHD_
    use flux_eos, only : source_b
#endif
    real(kind=DBL_KIND),intent(IN) :: dt
    integer :: n
    do n = MX, MX+NDIM-1
       call w_update_ndir(dt, n)
    end do
#ifdef _MHD_
    call source_b(F, W, dt)
#endif
  end subroutine w_update
  ! ---------------------------------------------------------------------------
  ! update v by flux for each direction
  ! ---------------------------------------------------------------------------
#define SHIFTR( A, NDIM ) cshift((A),  1, (NDIM) + DIMOFFSET)
#define SHIFTL( A, NDIM ) cshift((A), -1, (NDIM) + DIMOFFSET)
  subroutine w_update_ndir (dt, ndir)
    use grid
    real(kind=DBL_KIND),intent(IN) :: dt
    integer,intent(IN) :: ndir
    real(kind=DBL_KIND),dimension(MX:MZ) :: ds
    integer,parameter :: DIMOFFSET = 1-MX
    ds = get_ds()
    W = W - dt*ds(ndir)*(F(:,:,:,:,ndir) - SHIFTL(F(:,:,:,:,ndir),ndir))
  end subroutine w_update_ndir
  ! ---------------------------------------------------------------------------
  ! flux at cell interface for all directions
  ! ---------------------------------------------------------------------------
  subroutine get_flux
    integer :: n
    do n = MX, MX+NDIM-1
       call get_flux_ndir(n)
    end do
  end subroutine get_flux
  ! ---------------------------------------------------------------------------
  ! flux at cell interface for each direction
  ! ---------------------------------------------------------------------------
#define MINMOD(x, y) (max(0.d0,min((y)*sign(1.d0,(x)),abs(x)))*sign(1.d0,(x)))
#define SUPERBEE(x, y) (sign(1.d0,(y))*max(0.d0, min(sign(1.d0,(y))*BW*(x),abs(y)), min(sign(1.d0,(y))*(x),BW*abs(y))))
#if _MUSCL2_LIMITER_ == _MINMOD_
#define FLMT(x, y) MINMOD(x, y)
#endif
#if _MUSCL2_LIMITER_ == _SUPERBEE_
#define FLMT_(x, y) SUPERBEE(x, y)
#endif
  subroutine get_flux_ndir (ndir)
    use util
    use grid
    use flux_eos
    integer,intent(IN) :: ndir
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX) :: f1d, vl, vr
    integer,dimension(MMIN:MMAX) :: mcycle
    integer :: io,jo,ko,i2,j2,k2,i,j,k,m
#if _RECONSTRUCTION_ == _MUSCL2_
    real(kind=DBL_KIND),parameter :: BW = 2.d0
#endif
#if _RECONSTRUCTION_ == _MUSCL3_
    real(kind=DBL_KIND),parameter :: ETA = 1.d0/3.d0 ! 1/3 for 3rd order accuracy in space
    real(kind=DBL_KIND),parameter :: BW = (3.d0-ETA)/(1.d0-ETA)
!!$    real(kind=DBL_KIND),parameter :: ETA = -1.d0
!!$    real(kind=DBL_KIND),parameter :: BW = 1.d0
    real(kind=DBL_KIND) :: dva, dvb
#endif

    call util_arroffset(ndir,io,jo,ko)
    i2 = io*2
    j2 = jo*2
    k2 = ko*2
    do m = MMIN, MMAX
       do k = KMIN-ko, KMAX
          do j = JMIN-jo, JMAX
             do i = IMIN-io, IMAX
#if _RECONSTRUCTION_ == _NONE_
                vl(i,j,k,m) = V(i,j,k,m)
                vr(i,j,k,m) = V(i+io,j+jo,k+ko,m)
#endif
#if _RECONSTRUCTION_ == _MUSCL2_
                vl(i,j,k,m) = V(i,j,k,m) &
                     + (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i,j,k,m)-V(i-io,j-jo,k-ko,m)))*0.5d0
                vr(i,j,k,m) = V(i+io,j+jo,k+ko,m) &
                     - (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i+i2,j+j2,k+k2,m)-V(i+io,j+jo,k+ko,m)))*0.5d0
#endif
#if _RECONSTRUCTION_ == _MUSCL3_
                dva = V(i+io,j+jo,k+ko,m) - V(i,j,k,m)
                dvb = V(i,j,k,m) - V(i-io,j-jo,k-ko,m)
                vl(i,j,k,m) = V(i,j,k,m) &
                     + (1.d0 - ETA)/4.d0 * MINMOD(dvb, BW*dva) &
                     + (1.d0 + ETA)/4.d0 * MINMOD(dva, BW*dvb)
                dva = V(i+i2,j+j2,k+k2,m) - V(i+io,j+jo,k+ko,m)
                dvb = V(i+io,j+jo,k+ko,m) - V(i,j,k,m)
                vr(i,j,k,m) = V(i+io,j+jo,k+ko,m) &
                     - (1.d0 - ETA)/4.d0 * MINMOD(dva, BW*dvb) &
                     - (1.d0 + ETA)/4.d0 * MINMOD(dvb, BW*dva)
#endif
             enddo
          enddo
       enddo
    enddo
    mcycle = cyclecomp( ndir )
    vl = vl(:,:,:,mcycle)
    vr = vr(:,:,:,mcycle)
    call flux(vl, vr, f1d)
    F(:,:,:,mcycle,ndir) = f1d(:,:,:,:)

  end subroutine get_flux_ndir
end module timestep
