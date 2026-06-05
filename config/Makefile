.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# Paths
srcdir   = src
objdir   = include
libdir   = include
bindir   = bin
testdir  = tests

# Command-line options at make call
debug  ?= 0
openmp ?= 0

## COMPILER CONFIGURATION ##
# Assembled by configme from compiler -> machine -> netCDF fragments and
# spliced in at the placeholder below. The fragments define the compiler
# (FC, FFLAGS, FFLAGS_OPENMP, DFLAGS_*), netCDF wiring (INC_NC, LIB_NC),
# and any host-specific overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC, and LFLAGS_EXTRA set there. Defines the FESMUTILS
# and FFTW paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
    DFLAGS = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
    DFLAGS = $(DFLAGS_PROFILE)
endif

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
## Compilation rules (per-object; includes narrowed to what each file needs)
###############################################

$(objdir)/closures.o : $(srcdir)/closures.f90
	@mkdir -p $(objdir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c $< -o $@

$(objdir)/bucket.o : $(srcdir)/bucket.f90
	@mkdir -p $(objdir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c $< -o $@

$(objdir)/k24.o : $(srcdir)/k24.f90
	@mkdir -p $(objdir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_FFTW) -c $< -o $@

$(objdir)/fast_hydrology.o : $(srcdir)/fast_hydrology.f90 \
                             $(objdir)/closures.o \
                             $(objdir)/bucket.o \
                             $(objdir)/k24.o
	@mkdir -p $(objdir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) -c $< -o $@

###############################################
## Targets
###############################################

.PHONY: all fasthydro-static shmip greenland clean

all: fasthydro-static

fasthydro-static: $(fasthydro_objs)
	@mkdir -p $(libdir)
	ar rc $(libdir)/libfasthydro.a $(fasthydro_objs)
	ranlib $(libdir)/libfasthydro.a
	@echo ""
	@echo "    $(libdir)/libfasthydro.a is ready."
	@echo ""

shmip: fasthydro-static
	@mkdir -p $(bindir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) \
	    -o $(bindir)/shmip.x $(testdir)/shmip.f90 \
	    -L$(libdir) -lfasthydro $(LFLAGS)
	@echo ""
	@echo "    $(bindir)/shmip.x is ready."
	@echo ""

greenland: fasthydro-static
	@mkdir -p $(bindir)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) \
	    -o $(bindir)/greenland.x examples/greenland/greenland.f90 \
	    -L$(libdir) -lfasthydro $(LFLAGS)
	@echo ""
	@echo "    $(bindir)/greenland.x is ready."
	@echo ""

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(libdir)/libfasthydro.a $(bindir)/shmip.x $(bindir)/greenland.x
