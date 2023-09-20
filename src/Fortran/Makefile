#---------------------------------------------------
# OSX Monterey
#
MF91 = gfortran
MF90 = mpif90-openmpi-gcc12

FC1 = ${MF91} -I$(GACODE_ROOT)/modules -J$(GACODE_ROOT)/modules -fPIC -framework Accelerate
FC  = ${MF90} -I$(GACODE_ROOT)/modules -J$(GACODE_ROOT)/modules -fPIC -framework Accelerate -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
F77 = mpif77 -std=legacy
CC  = ${FC}

GOMPTEST := $(shell ${CC} -v 2>&1 | grep " version " | grep gcc >/dev/null; echo $$?)
ifeq ($(GOMPTEST),0)
  FOMP  =-fopenmp
endif

FMATH	=-fdefault-real-8 -fdefault-double-8 
FOPT	=-O3 
FDEBUG  =-Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -fbacktrace
F2PY    = f2py

# System math libraries

LMATH = /opt/local/lib/libfftw3.a
FFTW_INC = /opt/local/include 

ifdef FANN_ROOT
   # neural net libraries
   NN_LIB=-L$(GACODE_ROOT)/../neural/ -I$(GACODE_ROOT)/../neural/ -lbrainfuse -lfann
endif

# NetCDF libraries

NETCDF = -L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -lnetcdf -lnetcdf
NETCDF_INC = /opt/local/include
 
# Archive

ARCH = ar cr



########

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC1=tglf
EXEC2=tglf_mpi

LLIB=tglf_lib

ifdef FANN_ROOT
   # neural-net file
   NN=tglf_nn_TM.fann
else
   # stub file
   NN=tglf_nn_TM
endif

OBJECTS = tglf_isnan.o      \
          tglf_modules.o    \
          tglf_pkg.o        \
          tglf_allocate.o   \
          tglf_deallocate.o \
          tglf_startup.o    \
          tglf_hermite.o    \
          tglf_inout.o      \
          tglf_setup_geometry.o \
          tglf_LS.o         \
          tglf_eigensolver.o \
          tglf_geometry.o   \
          tglf_matrix.o     \
          tglf_max.o        \
          tglf_interface.o  \
          tglf_error.o      \
          tglf_isinf.o      \
          tglf_shutdown.o   \
          tglf_read_input.o \
          tglf_multiscale_spectrum.o \
          tglf_kygrid.o     \
          tglf_run.o        \
          tglf_TM.o         \
          tglf_init_mpi.o   \
          tglf_run_mpi.o    \
          tglf_TM_mpi.o     \
          $(NN).o

.SUFFIXES : .o .f90 .F90

all: $(LLIB).a $(EXEC1) Makefile

$(EXEC1): $(LLIB).a $(EXEC1).o
	$(FC1) $(FFLAGS) -o $(EXEC1) $(EXEC1).o $(LLIB).a $(LMATH) $(NN_LIB)

#$(EXEC1).so: $(LLIB).a $(EXEC1).o
#	$(FC1) $(FFLAGS) -shared -o $(EXEC1).so $(EXEC1).o $(LLIB).a $(LMATH) $(NN_LIB)

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

tglf_run_mpi.o: tglf_run.F90
	$(FC) -o tglf_run_mpi.o $(FMATH) $(FFLAGS) -DMPI_TGLF=1 -c $<

.f90.o:
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<
.F90.o:
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<

clean:
	rm -f *.so *.o *.lst $(LLIB).a $(EXEC1) $(EXEC2)
	rm -f ${GACODE_ROOT}/modules/tglf*.mod
