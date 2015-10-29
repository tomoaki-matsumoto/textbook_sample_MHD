module grid
  use parameter
  implicit none
  private
  ! V = primitive variables, U = conservative variables, W = conservative variables, F = flux
  real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX),target :: U, V, W
  real(kind=DBL_KIND),dimension(IMINGH:IMAXGH,JMINGH:JMAXGH,KMINGH:KMAXGH,MMIN:MMAX,MX:MX+NDIM-1),target :: F
  real(kind=DBL_KIND),dimension(IMINGH:IMAXGH),target :: X
  real(kind=DBL_KIND),dimension(JMINGH:JMAXGH),target :: Y
  real(kind=DBL_KIND),dimension(KMINGH:KMAXGH),target :: Z
  real(kind=DBL_KIND),dimension(MX:MZ) :: CellWidth
  real(kind=DBL_KIND) :: Time, Dtime
  integer :: Step
  logical,save :: Initialized = .FALSE.
  public :: U, V, W, F, X, Y, Z, Time, Dtime, Step, get_cellWidth, get_dv, get_ds
contains
  !-----------------------------------------------------------------------
  ! get cell width
  !-----------------------------------------------------------------------
  function get_cellWidth() result( r_cellwidth )
    real(kind=DBL_KIND),dimension(MX:MZ) :: r_cellwidth
    if ( .not. Initialized ) then
       CellWidth(MX) = X(IMIN+1)-X(IMIN)
       CellWidth(MY) = Y(JMIN+1)-Y(JMIN)
       CellWidth(MZ) = Z(KMIN+1)-Z(KMIN)
       Initialized = .TRUE.
    endif
    r_cellwidth(:) = CellWidth(:)
  end function get_cellWidth
  !-----------------------------------------------------------------------
  ! get dv
  !-----------------------------------------------------------------------
  function get_dv() result( dv )
    real(kind=DBL_KIND) :: dv
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    h = get_cellWidth()
    dv = h(MX)*h(MY)*h(MZ)
  end function get_dv
  !-----------------------------------------------------------------------
  ! get ds
  !-----------------------------------------------------------------------
  function get_ds() result( ds )
    real(kind=DBL_KIND),dimension(MX:MZ) :: ds
    real(kind=DBL_KIND),dimension(MX:MZ) :: h
    h = get_cellWidth()
    ds(MX) = h(MY)*h(MZ)
    ds(MY) = h(MX)*h(MZ)
    ds(MZ) = h(MX)*h(MY)
  end function get_ds

end module grid







