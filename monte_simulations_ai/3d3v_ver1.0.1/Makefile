# Makefile for 3D3V Non-Analog Track-Length Monte Carlo Simulation
#===============================================================================

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -fcheck=bounds -J$(MODDIR)
LDFLAGS =

# Faster production build (use: make FAST=1, optionally with OPENMP=1).
# Keep IEEE-sensitive transformations enabled by default; do not use -ffast-math.
ifdef FAST
  FFLAGS = -O3 -Wall -J$(MODDIR)
endif

# OpenMP support (use: make OPENMP=1)
ifdef OPENMP
  FFLAGS += -fopenmp
  LDFLAGS += -fopenmp
endif

# Directories
SRCDIR = code
RUNDIR = run
MODDIR = module
ifdef OPENMP
  MODDIR := $(MODDIR)_omp
endif
ifdef FAST
  MODDIR := $(MODDIR)_fast
endif

# Shared source files (in compilation order - dependencies first)
CORE_SOURCES = $(SRCDIR)/constants.f90 \
               $(SRCDIR)/data_types.f90 \
               $(SRCDIR)/random_utils.f90 \
               $(SRCDIR)/cdf_reader.f90 \
               $(SRCDIR)/tl_el_table_reader.f90 \
               $(SRCDIR)/cross_sections.f90 \
               $(SRCDIR)/source_terms.f90 \
               $(SRCDIR)/statistics.f90 \
               $(SRCDIR)/scoring.f90 \
               $(SRCDIR)/dynamics.f90 \
               $(SRCDIR)/io.f90

# Object files
CORE_OBJECTS = $(patsubst $(SRCDIR)/%.f90, $(MODDIR)/%.o, $(CORE_SOURCES))
MAIN_OBJECTS = $(CORE_OBJECTS) $(MODDIR)/main.o

# Executable
EXEC = $(RUNDIR)/monte_carlo_3d3v_natl

#===============================================================================
# Targets
#===============================================================================

.PHONY: all clean run plot plot-deltaE plot-energy plot-compare-deltaE force-link

all: $(EXEC)

$(EXEC): $(MAIN_OBJECTS) force-link
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(MAIN_OBJECTS) $(LDFLAGS)
	@echo "Build complete: $@"

force-link:

# Pattern rule for object files
$(MODDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(MODDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies
$(MODDIR)/data_types.o: $(MODDIR)/constants.o
$(MODDIR)/random_utils.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o
$(MODDIR)/cdf_reader.o: $(MODDIR)/constants.o
$(MODDIR)/tl_el_table_reader.o: $(MODDIR)/constants.o
$(MODDIR)/cross_sections.o: $(MODDIR)/constants.o $(MODDIR)/cdf_reader.o
$(MODDIR)/source_terms.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o
$(MODDIR)/statistics.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o
$(MODDIR)/scoring.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
                     $(MODDIR)/cross_sections.o $(MODDIR)/random_utils.o \
                     $(MODDIR)/cdf_reader.o $(MODDIR)/tl_el_table_reader.o \
                     $(MODDIR)/source_terms.o
$(MODDIR)/dynamics.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
                      $(MODDIR)/random_utils.o $(MODDIR)/cross_sections.o $(MODDIR)/cdf_reader.o
$(MODDIR)/io.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
                $(MODDIR)/random_utils.o $(MODDIR)/statistics.o
$(MODDIR)/main.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o $(MODDIR)/random_utils.o \
                  $(MODDIR)/cdf_reader.o $(MODDIR)/tl_el_table_reader.o \
                  $(MODDIR)/cross_sections.o $(MODDIR)/scoring.o \
                  $(MODDIR)/source_terms.o $(MODDIR)/statistics.o \
                  $(MODDIR)/dynamics.o $(MODDIR)/io.o

clean:
	rm -f module/*.o module/*.mod
	rm -f module_omp/*.o module_omp/*.mod
	rm -f module_fast/*.o module_fast/*.mod
	rm -f module_omp_fast/*.o module_omp_fast/*.mod
	rm -f $(EXEC)
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.mod
	@echo "Clean complete"

run: $(EXEC)
	cd $(RUNDIR) && ./monte_carlo_3d3v_natl

plot: plot-deltaE plot-energy

plot-deltaE:
	gnuplot plot_deltaE.gp

plot-energy:
	gnuplot plot_energy_hist.gp

plot-compare-deltaE:
	gnuplot plot_compare_deltaE.gp
