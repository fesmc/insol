# Shared build configuration for insol (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). This file references
# variables those provide: FFLAGS / FFLAGS_OPENMP (compiler) and LIB_NC
# (machine or auto-detected netCDF).

# fesm-utils: provides ncio, used by the irradiance test programs to read and
# write netCDF. The core library (libinsol.a) and its regression tests do not
# depend on it -- `make test` builds with nothing but a Fortran compiler.
FESMUTILSROOT = fesm-utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

# OpenMP build (make openmp=1): swap in the OpenMP variant of fesm-utils and
# append the compiler's OpenMP flag.
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

# Set by the machine fragment where the platform linker needs it; macOS ld
# rejects -zmuldefs, so the macbook fragment sets `LFLAGS_EXTRA =`.
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LIB_FESMUTILS) $(LIB_NC) $(LFLAGS_EXTRA)
