module io
  use parameter
  implicit none
  private
  integer,parameter :: STRLEN = 1000
  public :: io_writedata, STRLEN
contains
  subroutine io_writedata (dir)
    use grid
    integer,parameter :: LUN = 11
    character(LEN=*) :: dir
    character(LEN=STRLEN) :: filename
    call mkfilename( filename )
    print *, TRIM(filename)
    open(LUN, file=TRIM(dir)//'/'//TRIM(filename), form='unformatted')
    write(LUN) Time, Step
    write(LUN) IMAX, JMAX, KMAX, MMAX
    write(LUN) X, Y, Z
    write(LUN) V
    call flush(LUN)
    close(LUN)
  end subroutine io_writedata
  ! ---------------------------------------------------------------------------
  ! 
  ! ---------------------------------------------------------------------------
  subroutine mkfilename( filename )
    use grid
    character(LEN=STRLEN) :: filename
    character(LEN=2),parameter :: PREFIX = 'st'
    character(LEN=2),parameter :: SUFFIX = '.d'
    write(filename, '(A,I6.6,A)') PREFIX, Step, SUFFIX
  end subroutine mkfilename
end module io

