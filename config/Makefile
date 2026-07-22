.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir  = src
objdir  = libinsol/include
bindir  = libinsol/bin
testdir = tests

# Command-line options at make call
debug  ?= 0
openmp ?= 0

## COMPILER CONFIGURATION ##
# Assembled by `configme`: the compiler fragment (FC, FFLAGS, FFLAGS_OPENMP,
# DFLAGS_*), the machine fragment (host-specific overrides), and the detected
# netCDF wiring (INC_NC, LIB_NC).

<COMPILER_CONFIGURATION>

## DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS,
# FFLAGS_OPENMP, LIB_NC and LFLAGS_EXTRA set there.

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
##
## List of insol source files
##
###############################################

# Core library: orbital solution and insolation. No external dependencies.
insol_core = $(objdir)/insol_defs.o \
             $(objdir)/interp1D.o \
             $(objdir)/insolation.o

# Irradiance: surface irradiance decomposition (DNI/DHI) for the renewable
# energy applications. Independent of the core; see README.
insol_irradiance = $(objdir)/irradiance.o

insol_libs = $(insol_core) $(insol_irradiance)

$(objdir)/insol_defs.o: $(srcdir)/insol_defs.F90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(srcdir)/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(srcdir)/insolation.f90 $(objdir)/insol_defs.o $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/irradiance.o: $(srcdir)/irradiance.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
##
## Compilation of complete programs
##
###############################################

.PHONY: all static test test-dni check clean

all: static

# Static library
static: $(insol_libs)
	ar rc $(objdir)/libinsol.a $(insol_libs)
	ranlib $(objdir)/libinsol.a
	@echo " "
	@echo "    $(objdir)/libinsol.a is ready."
	@echo " "

# Regression tests for the core library. Deliberately free of netCDF and
# fesm-utils so the physics can be checked with only a Fortran compiler.
test: $(insol_core)
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_insol.x $^ $(testdir)/test_insol.f90
	@echo " "
	@echo "    $(bindir)/test_insol.x is ready."
	@echo " "

# Run the regression tests. Must be run from the repository root so that the
# input tables and reference values resolve.
check: test
	@$(bindir)/test_insol.x

# Irradiance / DNI tests. These read netCDF input, so they need fesm-utils.
test-dni: $(insol_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_NC) -o $(bindir)/test_dni.x \
		$^ $(testdir)/test_dni.f90 $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_dni.x is ready."
	@echo " "

clean:
	rm -f $(bindir)/*.x $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a
