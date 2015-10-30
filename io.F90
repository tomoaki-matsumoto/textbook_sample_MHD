#include "config.h"
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
    character(LEN=*),intent(IN) :: dir
#ifdef IO_TEXT
    call io_writedata_txt(dir)
#else
    call io_writedata_bin(dir)
#endif
  end subroutine io_writedata
  ! ---------------------------------------------------------------------------
  subroutine io_writedata_bin(dir)
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
  end subroutine io_writedata_bin
  ! ---------------------------------------------------------------------------
  subroutine io_writedata_txt(dir)
    use parameter
    character(LEN=*),intent(IN) :: dir
    if (NX /= 1 .and. NY == 1 .and. NZ ==1) then
       call io_writedata_txt_1d(dir)
    elseif (NX /= 1 .and. NY /= 1 .and. NZ ==1) then
       call io_writedata_txt_2d(dir)
    else
       print *, '*** This type of text output is not supported.'
    endif
  end subroutine io_writedata_txt
  ! ---------------------------------------------------------------------------
  subroutine io_writedata_txt_1d(dir)
    use grid
    integer,parameter :: LUN = 11
    character(LEN=*),intent(IN) :: dir
    character(LEN=STRLEN) :: filename
    integer :: i, m
    call mkfilename( filename, suffix='.txt')
    print *, 'Write '//TRIM(filename)
    open(LUN, file=TRIM(dir)//'/'//TRIM(filename))
    write(LUN, '("# time = ", E12.5)') Time
    write(LUN, '("# step = ", I6)') Step
    write(LUN, '("# imax = ", I6)') IMAX
    do i = 0, IMAX
       write(LUN, '(1PE13.5)', advance='no') X(i)
       do m = MMIN, MMAX
          write(LUN, '(1PE13.5)', advance='no') V(i,JMIN,KMIN,m)
       enddo
       write(LUN, '(A)') ''     ! newline
    enddo
    call flush(LUN)
    close(LUN)
  end subroutine io_writedata_txt_1d
  ! ---------------------------------------------------------------------------
  subroutine io_writedata_txt_2d(dir)
    use grid
    integer,parameter :: LUN = 11
    character(LEN=*),intent(IN) :: dir
    character(LEN=STRLEN) :: filename
    integer :: i, j, m
    call mkfilename( filename, suffix='.txt')
    print *, 'Write '//TRIM(filename)
    open(LUN, file=TRIM(dir)//'/'//TRIM(filename))
    write(LUN, '("# time = ", E12.5)') Time
    write(LUN, '("# step = ", I6)') Step
    write(LUN, '("# imax = ", I6)') IMAX
    write(LUN, '("# jmax = ", I6)') JMAX
    do j = 0, JMAX
       do i = 0, IMAX
          write(LUN, '(1PE13.5)', advance='no') X(i), Y(j)
          do m = MMIN, MMAX
             write(LUN, '(1PE13.5)', advance='no') V(i,j,KMIN,m)
          end do
          write(LUN, '(A)') ''     ! newline
       end do
    enddo
    call flush(LUN)
    close(LUN)
  end subroutine io_writedata_txt_2d
  ! ---------------------------------------------------------------------------
  subroutine mkfilename( filename, prefix, suffix )
    use grid
    character(LEN=STRLEN),intent(OUT) :: filename
    character(LEN=*),intent(IN),optional :: prefix, suffix
    character(LEN=2),parameter :: PREFIX_DEFAULT = 'st'
    character(LEN=2),parameter :: SUFFIX_DEFAULT = '.d'
    character(LEN=STRLEN) :: prefix_f, suffix_f
    if ( present(prefix) ) then
       prefix_f = prefix
    else
       prefix_f = PREFIX_DEFAULT
    endif
    if ( present(suffix) ) then
       suffix_f = suffix
    else
       suffix_f = SUFFIX_DEFAULT
    endif
    write(filename, '(A,I6.6,A)') trim(prefix_f), Step, trim(suffix_f)
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
