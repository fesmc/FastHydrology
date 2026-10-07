# Shared build configuration for FastHydrology (dependency wiring).
#
# Loaded after the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# Dependency paths (serial build by default).
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

FFTWROOT = fesm-utils/fftw/fftw-serial
INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib -lfftw3 -lm

# OpenMP build (make openmp=1): swap the serial deps for OpenMP variants and
# append the compiler's OpenMP flag (FFLAGS_OPENMP, set in the compiler fragment).
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    FFTWROOT = fesm-utils/fftw/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

    FFLAGS += $(FFLAGS_OPENMP)
endif

# Position-independent code (make pic=1): needed when the static library is
# linked into a shared library on Linux (e.g. yelmo's C API).
ifeq ($(pic), 1)
    FFLAGS += -fPIC
endif

LFLAGS_EXTRA ?=
LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_FFTW) $(LFLAGS_EXTRA)
