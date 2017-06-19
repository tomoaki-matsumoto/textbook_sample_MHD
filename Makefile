MAIN        = main
PARAMETER   = parameter
GRID        = grid
#BOUNDARY    = boundary_free
BOUNDARY    = boundary_periodic
#BOUNDARY    = boundary_windTunnel
#BOUNDARY    = boundary_bullet
#INIT        = init_wave
#INIT        = init_advect
#INIT        = init_shocktube
INIT        = init_Orszag_Tang
#INIT        = init_windTunnel
#INIT        = init_bullet
#FLUX        = flux_scalarAdvection
#FLUX        = flux_Roe
#FLUX        = flux_RoeM2
#FLUX        = flux_HLL
#FLUX        = flux_HLLC
FLUX        = flux_HLLD
#FLUX        = flux_HLLD_Boris
TIMESTEP    = timestep
UTIL        = util
IO          = io
CURE        = cure_crash
ERRNORM     = errornorm


#### intel fortran (ifort)
FC	 = ifort
FFLAGS = -u -O3 -shared-intel -mcmodel=large -fno-alias -fno-fnalias -openmp -openmp-report2
CPPFLAGS = 
# FFLAGS =  -traceback -g -warn all -check all -debug all

#### gnu gfotran 
# FC	 = gfortran
# FFLAGS = -O3 -ffree-line-length-none
# CPPFLAGS = 


OBJECT = \
	$(PARAMETER).o \
	$(UTIL).o \
	$(GRID).o \
	$(CURE).o \
	$(BOUNDARY).o \
	$(FLUX).o \
	$(TIMESTEP).o \
	$(INIT).o \
	$(IO).o \
	$(MAIN).o

.SUFFIXES:
.SUFFIXES: .o .f90 .F90
.f90.o:; $(FC) $(FFLAGS) $(CPPFLAGS) -c $<
.F90.o:; $(FC) $(FFLAGS) $(CPPFLAGS) -c $<

all : $(MAIN)

$(MAIN): $(OBJECT)
	$(FC) $(FFLAGS) -o $(MAIN) $(OBJECT)


OBJECT_ERRN = \
	$(PARAMETER).o \
	$(UTIL).o \
	$(GRID).o \
	$(IO).o \
	$(ERRNORM).o

$(ERRNORM): $(OBJECT_ERRN)
	$(FC) $(FFLAGS) -o $(ERRNORM) $(OBJECT_ERRN)

$(MAIN).o: $<

$(GRID).o: $<

$(PARAMETER).o: $< config.h

$(UTIL).o: $<

$(BOUNDARY).o: $< config.h

$(INIT).o: $< config.h

$(FLUX).o: $< config.h

$(TIMESTEP).o: $< config.h

$(IO).o: $< config.h

$(CURE).o: $< config.h

$(ERRNORM).o: $<

clean:
	-rm $(MAIN) $(OBJECT) $(ERRNORM) $(OBJECT_ERRN) *.mod

distclean:
	-rm $(MAIN) $(OBJECT) $(ERRNORM) $(OBJECT_ERRN) *.mod *.o

