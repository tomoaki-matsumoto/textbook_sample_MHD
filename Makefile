MAIN        = main
PARAMETER   = parameter
GRID        = grid
#BOUNDARY    = boundary_free
BOUNDARY    = boundary_periodic
#INIT        = init_advect
#INIT        = init_shocktube
INIT        = init_Orszag_Tang
#FLUX        = flux_Roe
#FLUX        = flux_HLLD
FLUX        = flux_HLLD_Boris
TIMESTEP    = timestep
UTIL        = util
IO          = io

FC	 = ifort
FFLAGS = -u -O3 -shared-intel -mcmodel=large -fno-alias -fno-fnalias
#FFLAGS =  -traceback -g -warn all -check all -debug all

OBJECT = \
	$(PARAMETER).o \
	$(UTIL).o \
	$(GRID).o \
	$(BOUNDARY).o \
	$(FLUX).o \
	$(TIMESTEP).o \
	$(INIT).o \
	$(IO).o \
	$(MAIN).o

.SUFFIXES:
.SUFFIXES: .o .f90 .F90
.f90.o:; $(FC) $(FFLAGS) -c $<
.F90.o:; $(FC) $(FFLAGS) -c $<

all : $(MAIN)

$(MAIN): $(OBJECT)
	$(FC) $(FFLAGS) -o $(MAIN) $(OBJECT)

$(MAIN).o: $<

$(GRID).o: $<

$(PARAMETER).o: $<

$(UTIL).o: $<

$(BOUNDARY).o: $<

$(INIT).o: $< config.h

$(FLUX).o: $<

$(TIMESTEP).o: $< config.h

$(IO).o: $<

clean:
	-rm $(MAIN) $(OBJECT) *.mod

distclean:
	-rm $(MAIN) $(OBJECT) *.mod *.o

