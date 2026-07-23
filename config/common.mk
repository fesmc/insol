# Shared build configuration for insol.
#
# Loaded after the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> common). insol has no external library
# dependencies -- no netCDF, no fesm-utils -- so this only handles the OpenMP
# build switch and the final link flags.

# OpenMP build (make openmp=1): append the compiler's OpenMP flag.
ifeq ($(openmp), 1)
    FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

# Set by the machine fragment where the platform linker needs it; macOS ld
# rejects -zmuldefs, so the macbook fragment sets `LFLAGS_EXTRA =`.
LFLAGS_EXTRA ?= -Wl,-zmuldefs

LFLAGS = $(LFLAGS_EXTRA)
