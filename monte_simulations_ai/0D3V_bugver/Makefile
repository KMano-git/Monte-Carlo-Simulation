# Makefile for 0D-3V Monte Carlo Simulation
#===============================================================================

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -fcheck=bounds -Jmodule
LDFLAGS =

# Directories
SRCDIR = code
RUNDIR = run

# Source files (in compilation order)
SOURCES = $(SRCDIR)/constants.f90 \
          $(SRCDIR)/particle_data.f90 \
          $(SRCDIR)/random_utils.f90 \
          $(SRCDIR)/cdf_reader.f90 \
          $(SRCDIR)/collisions.f90 \
          $(SRCDIR)/statistics.f90 \
          $(SRCDIR)/main.f90

# Object files
OBJECTS = $(patsubst $(SRCDIR)/%.f90, module/%.o, $(SOURCES))

# Module files
MODULES = constants.mod particle_data.mod random_utils.mod cdf_reader.mod \
          collisions.mod statistics.mod

# Execute
EXEC = $(RUNDIR)/monte_carlo_0d

#===============================================================================
# Targets
#===============================================================================

.PHONY: all clean run

all: $(EXEC)

$(EXEC): $(OBJECTS)
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Build complete: $@"

# Pattern rule for object files
module/%.o: $(SRCDIR)/%.f90
	@mkdir -p module
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies (module dependencies)
module/particle_data.o: module/constants.o
module/random_utils.o: module/constants.o
module/cdf_reader.o: module/constants.o
module/collisions.o: module/constants.o module/particle_data.o module/random_utils.o module/cdf_reader.o
module/statistics.o: module/constants.o module/particle_data.o
module/main.o: module/constants.o module/particle_data.o module/random_utils.o \
               module/cdf_reader.o module/collisions.o module/statistics.o

clean:
	rm -rf module/
	rm -f $(EXEC)
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.mod
	@echo "Clean complete"

run: $(EXEC)
	cd $(RUNDIR) && ./monte_carlo_0d
