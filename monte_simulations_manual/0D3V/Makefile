# Makefile for 0D-3V Monte Carlo Simulation
#===============================================================================

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -fcheck=bounds -J$(MODDIR)
LDFLAGS = 

# OpenMP support
OPENMP ?= 1
ifeq ($(OPENMP), 1)
    FFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif

# Directories
SRCDIR = code
RUNDIR = run
MODDIR = module

# Source files (in compilation order)
SOURCES = $(SRCDIR)/constants.f90 \
          $(SRCDIR)/particle_data.f90 \
          $(SRCDIR)/random_utils.f90 \
          $(SRCDIR)/cdf_reader.f90 \
		  $(SRCDIR)/cross_sections.f90 \
		  $(SRCDIR)/check_coll.f90 \
          $(SRCDIR)/collisions.f90 \
		  $(SRCDIR)/io.f90 \
          $(SRCDIR)/main.f90

# Object files
OBJECTS = $(patsubst $(SRCDIR)/%.f90, $(MODDIR)/%.o, $(SOURCES))

# Executable
EXEC = $(RUNDIR)/monte_carlo_0d3v

#===============================================================================
# Targets
#===============================================================================

.PHONY: all clean run

all: $(EXEC)

$(EXEC): $(OBJECTS)
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Build complete: $@"

# Pattern rule for object files (code/*.f90 → module/*.o)
$(MODDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies (module dependencies)
$(MODDIR)/particle_data.o: $(MODDIR)/constants.o
$(MODDIR)/random_utils.o: $(MODDIR)/constants.o $(MODDIR)/particle_data.o
$(MODDIR)/cdf_reader.o: $(MODDIR)/constants.o
$(MODDIR)/cross_sections.o: $(MODDIR)/constants.o $(MODDIR)/cdf_reader.o
$(MODDIR)/check_coll.o: $(MODDIR)/constants.o $(MODDIR)/particle_data.o $(MODDIR)/random_utils.o $(MODDIR)/cross_sections.o
$(MODDIR)/collisions.o: $(MODDIR)/constants.o $(MODDIR)/particle_data.o $(MODDIR)/random_utils.o $(MODDIR)/cdf_reader.o $(MODDIR)/cross_sections.o
$(MODDIR)/io.o: $(MODDIR)/constants.o $(MODDIR)/particle_data.o $(MODDIR)/random_utils.o
$(MODDIR)/main.o: $(MODDIR)/constants.o $(MODDIR)/particle_data.o $(MODDIR)/random_utils.o \
                  $(MODDIR)/cdf_reader.o $(MODDIR)/cross_sections.o $(MODDIR)/check_coll.o $(MODDIR)/collisions.o $(MODDIR)/io.o

clean:
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod
	rm -f $(EXEC)
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.mod
	@echo "Clean complete"

run: $(EXEC)
	cd $(RUNDIR) && ./monte_carlo_0d3v