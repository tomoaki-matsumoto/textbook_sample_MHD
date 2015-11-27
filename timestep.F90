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

#if defined(DIRECTIONAL_UNSPLIT) && defined(TIME_MARCHING_EULER)
    call step_unsplit_euler
#endif

#if defined(DIRECTIONAL_UNSPLIT) && defined(TIME_MARCHING_PC2)
    call step_unsplit_pc2
#endif

#if defined(DIRECTIONAL_UNSPLIT) && defined(TIME_MARCHING_RK2)
    call step_unsplit_rk2
#endif

#if defined(DIRECTIONAL_UNSPLIT) && defined(TIME_MARCHING_RK3)
    call step_unsplit_rk3
#endif

#if defined(DIRECTIONAL_UNSPLIT) && defined(TIME_MARCHING_K3)
    call step_unsplit_k3
#endif

#if defined(DIRECTIONAL_SPLIT)
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
    call boundary_fix(V)
  end subroutine step_unsplit_rk2
  ! ---------------------------------------------------------------------------
  ! 3rd order Runge-Kutta (Kutta method)
  ! Accuracy is same as that of RK3
  ! ---------------------------------------------------------------------------
  subroutine step_unsplit_k3
    use grid
    use flux_eos
    use boundary
    ! U^(n+1) = U^n + 1/6 (k1 + 4 k2 + k3)
    ! k1 = dt L^n
    ! k2 = dt L(U^n + k1 / 2)    = dt L(U*)
    ! k3 = dt L(U^n - k1 + 2 k2) = dt L(U**)
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),target :: k1, k2, k3
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux               ! k1
    W = 0
    call w_update(Dtime)        ! W := k1
    k1 = W
    W = U + k1/2                ! W := U*
    call u2v(W, V)              ! V := V*
    call boundary_fix(V)
    call get_flux               ! k2
    W = 0
    call w_update(Dtime)        ! W := k2
    k2 = W
    W = U - k1 + 2 * k2         ! W := U**
    call u2v(W, V)              ! V := V**
    call boundary_fix(V)
    call get_flux               ! k3
    W = 0
    call w_update(Dtime)        ! W := k3
    k3 = W
    W = U + (k1 + 4*k2 + k3)/6
    call u2v(W, V)
    call boundary_fix(V)
  end subroutine step_unsplit_k3
  ! ---------------------------------------------------------------------------
  ! 3rd order Runge-Kutta (see PLUTO users guide, Jiang & Shu (1996) JCP, 126(1), 202
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
    call boundary_fix(V)
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
    call boundary_fix(V)
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
    call boundary_fix(V)
  end subroutine step_unsplit_euler
  ! ---------------------------------------------------------------------------
  ! Macro for step_split_1d
  ! ---------------------------------------------------------------------------
#if defined(TIME_MARCHING_EULER)
#define STEP_SPLIT_1D(n) step_split_euler_1d(n)
#elif defined(TIME_MARCHING_PC2)
#define STEP_SPLIT_1D(n) step_split_pc2_1d(n)
#elif defined(TIME_MARCHING_RK2)
#define STEP_SPLIT_1D(n) step_split_rk2_1d(n)
#elif defined(TIME_MARCHING_RK3)
#define STEP_SPLIT_1D(n) step_split_rk3_1d(n)
#elif defined(TIME_MARCHING_K3)
#define STEP_SPLIT_1D(n) step_split_k3_1d(n)
#else
  ERROR
#endif
  ! ---------------------------------------------------------------------------
  ! fractional time step
  ! ---------------------------------------------------------------------------
  subroutine step_split
#ifdef FLUX_SCHEME_MHD
    use grid, only: V, W, F, Dtime
    use flux_eos, only : source_b, v2u, u2v
#endif !FLUX_SCHEME_MHD
    if (NDIM == 1) then
       call STEP_SPLIT_1D(MX)
    elseif (NDIM == 2) then
       call step_split_2d
    elseif (NDIM == 3) then
       call step_split_3d
    else
       print *, '*** error in step_fractional'
    endif
#ifdef FLUX_SCHEME_MHD
    call v2u(V, W)
    call source_b(F, W, Dtime)
    call u2v(W, V)
#endif !FLUX_SCHEME_MHD
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
    call boundary_fix(V)
  end subroutine step_split_rk2_1d
  ! ---------------------------------------------------------------------------
  ! 3rd order Runge-Kutta (Kutta method)
  ! ---------------------------------------------------------------------------
  subroutine step_split_k3_1d(ndir)
    use grid
    use flux_eos
    use boundary
    integer,intent(IN) :: ndir
    ! U^(n+1) = U^n + 1/6 (k1 + 4 k2 + k3)
    ! k1 = dt L^n
    ! k2 = dt L(U^n + k1 / 2)    = dt L(U*)
    ! k3 = dt L(U^n - k1 + 2 k2) = dt L(U**)
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),target :: k1, k2, k3
    call v2u(V, U)
    call boundary_fix(V)
    call get_flux_ndir(ndir)    ! k1
    W = 0
    call w_update_ndir(Dtime,ndir) ! W := k1
    k1 = W
    W = U + k1/2                ! W := U*
    call u2v(W, V)              ! V := V*
    call boundary_fix(V)
    call get_flux_ndir(ndir)    ! k2
    W = 0
    call w_update_ndir(Dtime,ndir) ! W := k2
    k2 = W
    W = U - k1 + 2 * k2         ! W := U**
    call u2v(W, V)              ! V := V**
    call boundary_fix(V)
    call get_flux_ndir(ndir)    ! k3
    W = 0
    call w_update_ndir(Dtime,ndir) ! W := k3
    k3 = W
    W = U + (k1 + 4*k2 + k3)/6
    call u2v(W, V)
    call boundary_fix(V)
  end subroutine step_split_k3_1d
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
    call boundary_fix(V)
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
    call boundary_fix(V)
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
    call boundary_fix(V)
  end subroutine step_split_euler_1d
  ! ---------------------------------------------------------------------------
  ! update u by flux for all directions
  ! ---------------------------------------------------------------------------
  subroutine w_update( dt )
    use grid
#ifdef FLUX_SCHEME_MHD
    use flux_eos, only : source_b
#endif !FLUX_SCHEME_MHD
    real(kind=DBL_KIND),intent(IN) :: dt
    integer :: n
    do n = MX, MX+NDIM-1
       call w_update_ndir(dt, n)
    end do
#ifdef FLUX_SCHEME_MHD
    call source_b(F, W, dt)
#endif !FLUX_SCHEME_MHD
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
#ifdef MUSCL2_LIMITER_MINMOD
#define FLMT(x, y) MINMOD(x, y)
#endif
#ifdef MUSCL2_LIMITER_SUPERBEE
#define FLMT(x, y) SUPERBEE(x, y)
#endif
#ifdef MUSCL2_WO_LIMITER
#define FLMT(x, y) (y)
#endif
  subroutine get_flux_ndir (ndir)
    use util
    use grid
    use flux_eos
    integer,intent(IN) :: ndir
    real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX) :: f1d, vl, vr
    integer,dimension(MMIN:MMAX) :: mcycle
    integer :: io,jo,ko,i2,j2,k2,i,j,k,m
#ifdef RECONSTRUCTION_MUSCL2
    real(kind=DBL_KIND),parameter :: BW = 2.d0
#endif
#ifdef RECONSTRUCTION_MUSCL3
    real(kind=DBL_KIND),parameter :: ETA = 1.d0/3.d0 ! 1/3 for 3rd order accuracy in space
    real(kind=DBL_KIND),parameter :: BW = (3.d0-ETA)/(1.d0-ETA)
!!$    real(kind=DBL_KIND),parameter :: ETA = -1.d0
!!$    real(kind=DBL_KIND),parameter :: BW = 1.d0
    real(kind=DBL_KIND) :: dva, dvb
#endif
#ifdef RECONSTRUCTION_LIMO3
    real(kind=DBL_KIND) :: thtL, thtR, duLL, duLR, duRL, duRR, phiL, phiR, etaL, etaR, flagL, flagR
    real(kind=DBL_KIND),parameter :: rlimit = 1.d0
    real(kind=DBL_KIND),parameter :: eps = 1.d-20
    real(kind=DBL_KIND) :: cellWidth(MX:MZ), dx
#endif

    call util_arroffset(ndir,io,jo,ko)
    i2 = io*2
    j2 = jo*2
    k2 = ko*2
#ifdef RECONSTRUCTION_LIMO3
    cellWidth=get_cellWidth()
    dx = cellWidth(ndir)
#endif
    do m = MMIN, MMAX
       do k = KMIN-ko, KMAX
          do j = JMIN-jo, JMAX
             do i = IMIN-io, IMAX
#if defined(RECONSTRUCTION_NONE)
                vl(i,j,k,m) = V(i,j,k,m)
                vr(i,j,k,m) = V(i+io,j+jo,k+ko,m)
#elif defined(RECONSTRUCTION_MUSCL2)
                vl(i,j,k,m) = V(i,j,k,m) &
                     + (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i,j,k,m)-V(i-io,j-jo,k-ko,m)))*0.5d0
                vr(i,j,k,m) = V(i+io,j+jo,k+ko,m) &
                     - (FLMT(V(i+io,j+jo,k+ko,m)-V(i,j,k,m), V(i+i2,j+j2,k+k2,m)-V(i+io,j+jo,k+ko,m)))*0.5d0
#elif defined(RECONSTRUCTION_MUSCL3)
                ! Computational Gasdynamics, p 581, C. B. Laney (1998)
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
#elif defined(RECONSTRUCTION_LIMO3)
                ! Cada & Torrilhon (2009), JCP, 228, 4118, 10.1016/j.jcp.2009.02.020
                duLL = V(i,j,k,m)-V(i-io,j-jo,k-ko,m)
                duLR = V(i+io,j+jo,k+ko,m)-V(i,j,k,m)
                duRL = duLR
                duRR = V(i+i2,j+j2,k+k2,m)-V(i+io,j+jo,k+ko,m)
                thtL = duLL/(abs(duLR)+eps)*sign(1.d0,duLR)
                thtR = duRR/(abs(duRL)+eps)*sign(1.d0,duRL)
                etaL = (duLL**2 + duLR**2)/(rlimit*dx)**2
                etaR = (duRL**2 + duRR**2)/(rlimit*dx)**2
                flagL = 0.5d0+sign(0.5d0,etaL-1.d0) ! if eta > 1 then 1 else 0
                flagR = 0.5d0+sign(0.5d0,etaL-1.d0) ! if eta > 1 then 1 else 0
                phiL = flagL * max(0.d0, min((2.d0+thtL)/3.d0, max(-0.5d0*thtL,min(2.d0*thtL,(2.d0+thtL)/3.d0,1.6d0)))) &
                     + (1.d0-flagL) * (2.d0+thtL)/3.d0
                phiR = flagR * max(0.d0, min((2.d0+thtR)/3.d0, max(-0.5d0*thtR,min(2.d0*thtR,(2.d0+thtR)/3.d0,1.6d0)))) &
                     + (1.d0-flagR) * (2.d0+thtR)/3.d0
                vl(i,j,k,m) = V(i,j,k,m) + 0.5d0*duLR*phiL
                vr(i,j,k,m) = V(i+io,j+jo,k+ko,m) - 0.5d0*duRL*phiR
#else
                ERROR
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
