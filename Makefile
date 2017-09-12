OUTPUT_DIR ?= build
USETMP ?= 1

BUILD ?= MPI

#dimension of simulation space; should be either 2 or 3
#DIM=2
DIM=3

#velocity estimation; 0 for baseline; 1 for Yamaguchi & Yokoyama 2002
VEL=0

ifortErr = $(shell which ifort >/dev/null; echo $$?)
#these settings for ifort 13 and higher
ifeq "$(ifortErr)" "0"
#ifort; Can remove -xHost if your cluster is not uniform, or specify specific processor optimizations -x...
F90C     = ifort
MPIF90C = mpiifort
### optimization?
#!#!#FFLAGS = -mt_mpi -parallel -qopenmp -fpp -mkl:cluster -qopt-report5 -xHOST -ipo -no-prec-div -fp-model fast=2 -static-intel -align array64byte
#FFLAGS = -mt_mpi -qopenmp -fpp -mkl:cluster -qopt-report5 -xHOST -O3 -ipo -no-prec-div -fp-model fast=2 -static-intel
#FFLAGS = -mt_mpi -qopenmp -fpp -mkl:cluster -qopt-report5 -xHOST -O3 -ipo -no-prec-div -fp-model fast=2 -static-intel -i8
### debugging
FFLAGS = -mt_mpi -mkl:cluster -O0 -static-intel -qopenmp -g -check all -check noarg_temp_created -traceback -no-prec-div -fpp -fpe0 -debug
#FFLAGS = -mt_mpi -qopenmp -fpp -mkl:cluster -qopt-report5 -xHOST -O0 -ipo -no-prec-div -fp-model fast=2 -static-intel -g -check all -check noarg_temp_created -traceback -fpe0 -debug
#FFLAGS = -mt_mpi -qopenmp -fpp -mkl:cluster -qopt-report5 -xHOST -O0 -ipo -no-prec-div -fp-model fast=2 -static-intel -i8 -static-intel -g -check all -check noarg_temp_created -traceback -fpe0 -debug
MODOUT = -module $(OUTPUT_DIR)
LAPACKL =
#FFTWL = -L$(FFTW_ROOT)/lib -lfftw3_mpi -lfftw3 -lm
#FFTWL = -L$(MKLROOT)/lib -lfftw3x_cdft_ilp64
#FFTWL = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl -L./fftw3x_cdft -lfftw3x_cdft_ilp64
#FFTWL = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -lm -ldl -L./fftw3x_cdft -lfftw3x_cdft_ilp64
FFTWL = -L./fftw3x_cdft -lfftw3x_cdft_ilp64

else
endif

#use PRECISION=SINGLE to use single precision (probably doesn't work)
PRECISION ?= 

#INCLUDE = -I$(FFTW_ROOT)/include
INCLUDE = -I$(MKLROOT)/include/fftw
#?#INCLUDE = -I$(MKLROOT)/include

F90FLAGS = $(FFLAGS) -DDIM=$(DIM) $(INCLUDE) -DVEL=$(VEL) -DUSETMP=$(USETMP)
#!#LINKFLAGS = $(LAPACKL) $(FFTWL)
LINKFLAGS = $(LAPACKL) $(FFTWL) -liomp5 -lpthread

ifeq ($(VEL),1)
SUFFIX_VEL=_velYY
else
SUFFIX_VEL=
endif

ifneq ($(PRECISION),)
FFLAGS += -D$(PRECISION) -DMATRIX_$(PRECISION)
endif

ifeq ($(BUILD),MPI)
override OUTPUT_DIR :=$(OUTPUT_DIR)_MPI
FFLAGS +=  -DMPI
LINKFLAGS +=  $(LINKMPI)
F90C = $(MPIF90C)
endif

ifeq ($(OFFLOAD),mandatory)
SUFFIX_OFFLOAD = _offload
override OUTPUT_DIR :=$(OUTPUT_DIR)$(SUFFIX_OFFLOAD)
endif

# routines related to defects
ifeq "$(DIM)" "2"
#DEFECTS=defects_2d.o
DEFECTS=defects_none.o
else ifeq "$(DIM)" "3"
DEFECTS=defects_3d.o
else
# error message for DIM/=2,3???
endif

OBJFILES = $(OUTPUT_DIR)/utils.o $(OUTPUT_DIR)/subroutines.o $(OUTPUT_DIR)/inifile.o $(OUTPUT_DIR)/potential.o $(OUTPUT_DIR)/settings.o $(OUTPUT_DIR)/$(DEFECTS) $(OUTPUT_DIR)/initial.o $(OUTPUT_DIR)/spectrum.o $(OUTPUT_DIR)/grav_waves.o $(OUTPUT_DIR)/evolve.o $(OUTPUT_DIR)/main.o

default: axion

$(OUTPUT_DIR)/subroutines.o: $(OUTPUT_DIR)/utils.o
$(OUTPUT_DIR)/settings.o: $(OUTPUT_DIR)/subroutines.o $(OUTPUT_DIR)/potential.o
$(OUTPUT_DIR)/$(DEFECTS): $(OUTPUT_DIR)/settings.o
$(OUTPUT_DIR)/initial.o: $(OUTPUT_DIR)/$(DEFECTS)
$(OUTPUT_DIR)/spectrum.o: $(OUTPUT_DIR)/initial.o
$(OUTPUT_DIR)/grav_wavs.o: $(OUTPUT_DIR)/spectrum.o
$(OUTPUT_DIR)/evolve.o: $(OUTPUT_DIR)/grav_waves.o
$(OUTPUT_DIR)/main.o: $(OUTPUT_DIR)/evolve.o

F90FLAGS += $(MODOUT) -I$(OUTPUT_DIR)/

#export FFLAGS
#export F90C
#export OUTPUT_DIR

directories:
	mkdir -p $(OUTPUT_DIR)

$(OUTPUT_DIR)/%.o: %.c
	$(CC) $(GSLINC) -c $*.c -o $(OUTPUT_DIR)/$*.o

$(OUTPUT_DIR)/%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90 -o $(OUTPUT_DIR)/$*.o

axion: directories $(OBJFILES)
	$(F90C) -o ../axion_$(DIM)d$(SUFFIX_VEL)$(SUFFIX_OFFLOAD) $(OBJFILES) $(LINKFLAGS) $(F90FLAGS)

clean:
	rm -f $(OUTPUT_DIR)/* ../core
