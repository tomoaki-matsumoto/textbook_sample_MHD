!
! print L1norm, L2norm, Linfnorm, evaluating data in DATA directory
!
program main
  use parameter
  use grid
  use io
  implicit none
  character(LEN=STRLEN) :: dir, path
  real(kind=DBL_KIND),dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MMIN:MMAX) :: diff, q0, q1
  real(kind=DBL_KIND) :: normL1, normL2, normLinf
  integer :: i
  dir = 'DATA'
  call io_unix_command('ls '//TRIM(dir)//' |head -1', path)
  call io_readdata(TRIM(dir)//'/'//TRIM(path))
  q0 = V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MMIN:MMAX)
  call io_unix_command('ls '//TRIM(dir)//' |tail -1', path)
  call io_readdata(TRIM(dir)//'/'//TRIM(path))
  q1 = V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,MMIN:MMAX)
  diff = q1 - q0

  normL1 = sum(abs(diff))/size(diff)
  normL2 = sqrt(sum(diff**2)/size(diff))
  normLinf = maxval(diff)

  write(6, '(I5,3E12.5)') NX, normL1, normL2, normLinf

end program main
