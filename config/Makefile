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
# DFLAGS_*) and the machine fragment (host-specific overrides). insol has no
# external library dependencies, so no netCDF wiring is needed.

<COMPILER_CONFIGURATION>

## DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference FFLAGS and
# FFLAGS_OPENMP set there.

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
             $(objdir)/insol_orbit.o \
             $(objdir)/insolation.o

# Irradiance: surface irradiance decomposition (DNI/DHI) for the renewable
# energy applications. Shares the orbital solution (insol_orbit) with the core.
insol_irradiance = $(objdir)/irradiance.o

insol_libs = $(insol_core) $(insol_irradiance)

$(objdir)/insol_defs.o: $(srcdir)/insol_defs.F90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insol_orbit.o: $(srcdir)/insol_orbit.f90 $(objdir)/insol_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(srcdir)/insolation.f90 $(objdir)/insol_orbit.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/irradiance.o: $(srcdir)/irradiance.f90 $(objdir)/insol_orbit.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
##
## Compilation of complete programs
##
###############################################

.PHONY: all static test check clean

all: static

# Static library
static: $(insol_libs)
	ar rc $(objdir)/libinsol.a $(insol_libs)
	ranlib $(objdir)/libinsol.a
	@echo " "
	@echo "    $(objdir)/libinsol.a is ready."
	@echo " "

# Regression tests. The whole library builds with only a Fortran compiler.
test: $(insol_core) $(insol_irradiance)
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_insol.x $(insol_core) $(testdir)/test_insol.f90
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_orbit.x $(insol_core) $(testdir)/test_orbit.f90
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_irradiance.x $(objdir)/insol_defs.o $(objdir)/insol_orbit.o $(insol_irradiance) $(testdir)/test_irradiance.f90
	@echo " "
	@echo "    test_insol.x, test_orbit.x, test_irradiance.x are ready."
	@echo " "

# Run the regression tests. Must be run from the repository root so that the
# input tables and reference values resolve.
check: test
	@$(bindir)/test_insol.x
	@$(bindir)/test_orbit.x
	@$(bindir)/test_irradiance.x

clean:
	rm -f $(bindir)/*.x $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a
