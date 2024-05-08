#-----------------------------------------------------
# Makefile to compile the TGLFEP_driver system.
#-----------------------------------------------------

# Define compilers and flags.
# Compilers and flags

include ${GACODE_ROOT}/platform/build/make.inc.${GACODE_PLATFORM}
export EXTRA_LIBS = \
        ${GACODE_ROOT}/tglf/src/tglf_lib.a \
        ${GACODE_ROOT}/f2py/expro/expro_lib.a \
        ${GACODE_ROOT}/f2py/geo/geo_lib.a \
        ${GACODE_ROOT}/shared/math/math_lib.a \
        ${GACODE_ROOT}/shared/UMFPACK/UMFPACK_lib.a
ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=TGLFEP_driver

LLIB=TGLFEP_lib

OBJECTS = TGLFEP_interface.o\
          TGLFEP_tglf_map.o\
          TGLFEP_read_EXPRO.o\
          TGLFEP_ky.o\
          TGLFEP_TM.o\
          TGLFEP_ky_widthscan.o\
          TGLFEP_ky_nEPscan.o\
          TGLFEP_scalefactor.o\
	  TGLFEP_kwscale_scan.o\
          TGLFEP_interpolate.o\
          TGLFEP_complete_output.o\
          TGLFEP_EPflux.o\
          TGLFEP_mainsub.o 

.SUFFIXES : .o .f90 .f

all: $(LLIB).a $(EXEC)

$(EXEC): $(LLIB).a $(EXEC).o ${EXTRA_LIBS}
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o $(LLIB).a $(EXTRA_LIBS) ${LMATH} $(NN_LIB)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) $(NN_LIB) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXEC) $(LLIB).a
