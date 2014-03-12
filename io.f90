module io
  use parameter
  implicit none
  private
  integer,parameter :: STRLEN = 1000
  public :: io_writedata, io_readdata, io_unix_command, STRLEN
contains
  ! ---------------------------------------------------------------------------
  subroutine io_readdata(fullpath)
    use grid
    character(LEN=*),intent(IN) :: fullpath
    integer,parameter :: LUN = 11
    integer :: imax_r, jmax_r, kmax_r, mmax_r
    print *, 'Read '//TRIM(fullpath)
    open(LUN, file=TRIM(fullpath), form='unformatted')
    read(LUN) Time, Step
    read(LUN) imax_r, jmax_r, kmax_r, mmax_r
    read(LUN) X, Y, Z
    read(LUN) V
    close(LUN)
    if (imax_r /= IMAX) then 
       print *, '** error in IMAX'
       stop
    endif
    if (jmax_r /= JMAX) then 
       print *, '** error in JMAX'
       stop
    endif
    if (kmax_r /= KMAX) then 
       print *, '** error in KMAX'
       stop
    endif
    if (mmax_r /= MMAX) then 
       print *, '** error in MMAX'
       stop
    endif
  end subroutine io_readdata
  ! ---------------------------------------------------------------------------
  subroutine io_writedata(dir)
    use grid
    integer,parameter :: LUN = 11
    character(LEN=*),intent(IN) :: dir
    character(LEN=STRLEN) :: filename
    call mkfilename( filename )
    print *, 'Write '//TRIM(filename)
    open(LUN, file=TRIM(dir)//'/'//TRIM(filename), form='unformatted')
    write(LUN) Time, Step
    write(LUN) IMAX, JMAX, KMAX, MMAX
    write(LUN) X, Y, Z
    write(LUN) V
    call flush(LUN)
    close(LUN)
  end subroutine io_writedata
  ! ---------------------------------------------------------------------------
  subroutine mkfilename( filename )
    use grid
    character(LEN=STRLEN),intent(OUT) :: filename
    character(LEN=2),parameter :: PREFIX = 'st'
    character(LEN=2),parameter :: SUFFIX = '.d'
    write(filename, '(A,I6.6,A)') PREFIX, Step, SUFFIX
  end subroutine mkfilename
  ! ---------------------------------------------------------------------------
  subroutine io_unix_command(command, buf)
    character(LEN=*),intent(IN) :: command
    character(LEN=*),intent(OUT) :: buf
    character(LEN=STRLEN),parameter :: TMPFILE = '/tmp/io_unix_command'
    integer :: i, system
    integer,parameter :: LUN=11
    i = system( TRIM(command)//">"//TRIM(TMPFILE)//" 2>&1" )
    open(LUN, file=TMPFILE)
    read(LUN, *) buf
    call flush(LUN)
    close(LUN)
    i = system( "/bin/rm -f "//TRIM(TMPFILE) )
  end subroutine io_unix_command
end module io
