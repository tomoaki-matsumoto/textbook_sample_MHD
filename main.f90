program main
  use parameter
  use grid
  use io
  use init
  use flux_eos
  use timestep
  implicit none
  integer,parameter :: STEPMAX = 1000000
!!$  real(kind=DBL_KIND),parameter :: T_LAST = 0.2d0
  real(kind=DBL_KIND),parameter :: T_LAST = PI
!!$  real(kind=DBL_KIND),parameter :: T_LAST = 1.d0
  logical :: bool_halt = .FALSE.
  character(LEN=STRLEN) :: dir
  Step=0
  dir = 'DATA'
  call system("rm -rf "//TRIM(dir))
  call system("mkdir -p "//TRIM(dir))
  call init_cond
  call io_writedata(dir)
  do 
     call cflcond
     if (Time + Dtime > T_LAST) then
        Dtime = T_LAST - Time
        bool_halt = .true.
     endif
     Time = Time + Dtime
     Step = Step + 1
     print *, Step, Dtime, Time
     call step_full
     if ( mod(Step, 10) == 0 ) call io_writedata(dir)
     if ( Step >= STEPMAX ) exit
     if ( bool_halt ) exit
  end do
  call io_writedata(dir)
end program main




