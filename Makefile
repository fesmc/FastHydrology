.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Paths
srcdir   = src
objdir   = include
libdir   = include
bindir   = bin
testdir  = tests

# Build options
debug  ?= 0
openmp ?= 0

# Compiler
FC = gfortran

FFLAGS_BASE   = -ffree-line-length-none -I$(objdir) -J$(objdir)
FFLAGS        = $(FFLAGS_BASE)
FFLAGS_OPENMP = -fopenmp

DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all

DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
    DFLAGS = $(DFLAGS_DEBUG)
endif

# fesm-utils (provides nml, FFTW, LIS). Override FESMUTILS_DIR on the make
# command line if your checkout lives elsewhere, e.g.
#   make FESMUTILS_DIR=/path/to/fesm-utils
FESMUTILS_DIR ?= ../yelmo/fesm-utils
FESMUTILSROOT  = $(FESMUTILS_DIR)/utils

ifeq ($(openmp), 1)
    INC_FESMUTILS = -I$(FESMUTILSROOT)/include-omp
    LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-omp -lfesmutils
    FFTWROOT      = $(FESMUTILS_DIR)/fftw-omp
    INC_FFTW      = -I$(FFTWROOT)/include
    LIB_FFTW      = -L$(FFTWROOT)/lib -lfftw3_omp -lfftw3 -lm
    FFLAGS += $(FFLAGS_OPENMP)
else
    INC_FESMUTILS = -I$(FESMUTILSROOT)/include-serial
    LIB_FESMUTILS = -L$(FESMUTILSROOT)/include-serial -lfesmutils
    FFTWROOT      = $(FESMUTILS_DIR)/fftw-serial
    INC_FFTW      = -I$(FFTWROOT)/include
    LIB_FFTW      = -L$(FFTWROOT)/lib -lfftw3 -lm
endif

# NetCDF (autodetected via nf-config / nc-config if available)
INC_NC ?= $(shell nf-config --fflags 2>/dev/null) $(shell nc-config --cflags 2>/dev/null)
LIB_NC ?= $(shell nf-config --flibs  2>/dev/null) $(shell nc-config --libs   2>/dev/null)

LFLAGS_EXTRA ?=
LFLAGS = $(LIB_FESMUTILS) $(LIB_FFTW) $(LIB_NC) $(LFLAGS_EXTRA)

###############################################
## Source list
###############################################

fasthydro_src = \
    $(srcdir)/closures.f90 \
    $(srcdir)/bucket.f90 \
    $(srcdir)/k24.f90 \
    $(srcdir)/fast_hydrology.f90

fasthydro_objs = \
    $(objdir)/closures.o \
    $(objdir)/bucket.o \
    $(objdir)/k24.o \
    $(objdir)/fast_hydrology.o

###############################################
## Compilation rules
###############################################

$(objdir)/%.o : $(srcdir)/%.f90
	@mkdir -p $(objdir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_FFTW) $(INC_NC) -c $< -o $@

# Explicit dependencies (Fortran module order)
$(objdir)/closures.o       : $(srcdir)/closures.f90
$(objdir)/bucket.o         : $(srcdir)/bucket.f90
$(objdir)/k24.o            : $(srcdir)/k24.f90
$(objdir)/fast_hydrology.o : $(srcdir)/fast_hydrology.f90 \
                             $(objdir)/closures.o \
                             $(objdir)/bucket.o \
                             $(objdir)/k24.o

###############################################
## Targets
###############################################

.PHONY: all lib shmip greenland clean

all: lib

lib: $(fasthydro_objs)
	@mkdir -p $(libdir)
	ar rc $(libdir)/libfasthydro.a $(fasthydro_objs)
	ranlib $(libdir)/libfasthydro.a
	@echo ""
	@echo "    $(libdir)/libfasthydro.a is ready."
	@echo ""

shmip: lib
	@mkdir -p $(bindir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_FFTW) $(INC_NC) \
	    -o $(bindir)/shmip.x $(testdir)/shmip.f90 \
	    -L$(libdir) -lfasthydro $(LFLAGS)
	@echo ""
	@echo "    $(bindir)/shmip.x is ready."
	@echo ""

greenland: lib
	@mkdir -p $(bindir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_FFTW) $(INC_NC) \
	    -o $(bindir)/greenland.x examples/greenland/greenland.f90 \
	    -L$(libdir) -lfasthydro $(LFLAGS)
	@echo ""
	@echo "    $(bindir)/greenland.x is ready."
	@echo ""

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(libdir)/libfasthydro.a $(bindir)/shmip.x $(bindir)/greenland.x
