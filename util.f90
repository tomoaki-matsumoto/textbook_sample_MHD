module util
  implicit none
contains
  !-----------------------------------------------------------------------
  ! arrey offset
  ! when ncrd = MX, io = 1
  ! when ncrd = MY, jo = 1
  ! when ncrd = MZ, jo = 1
  !-----------------------------------------------------------------------
  subroutine util_arroffset(ncrd,io,jo,ko)
    use parameter
    integer,intent(IN) :: ncrd
    integer,intent(OUT) :: io,jo,ko
    io = 0
    jo = 0
    ko = 0
    select case(ncrd)
    case ( MX )                 ! x-direction
       io = 1
    case ( MY )                 ! y-direction
       jo = 1
    case ( MZ )                 ! z-direction
       ko = 1
    case default
       write(*,*) '*** bad ncrd in subroutine arroffset'
       stop
    end select
  end subroutine util_arroffset
end module util
