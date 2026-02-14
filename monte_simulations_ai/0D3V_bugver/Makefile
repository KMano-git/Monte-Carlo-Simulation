# Makefile for 0D-3V Monte Carlo Simulation
#===============================================================================

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -fcheck=bounds
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
OBJECTS = $(SOURCES:.f90=.o)

# Module files
MODULES = constants.mod particle_data.mod random_utils.mod cdf_reader.mod \
          collisions.mod statistics.mod

# Executable
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
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies (module dependencies)
$(SRCDIR)/particle_data.o: $(SRCDIR)/constants.o
$(SRCDIR)/random_utils.o: $(SRCDIR)/constants.o
$(SRCDIR)/cdf_reader.o: $(SRCDIR)/constants.o
$(SRCDIR)/collisions.o: $(SRCDIR)/constants.o $(SRCDIR)/particle_data.o $(SRCDIR)/random_utils.o $(SRCDIR)/cdf_reader.o
$(SRCDIR)/statistics.o: $(SRCDIR)/constants.o $(SRCDIR)/particle_data.o
$(SRCDIR)/main.o: $(SRCDIR)/constants.o $(SRCDIR)/particle_data.o $(SRCDIR)/random_utils.o \
                  $(SRCDIR)/cdf_reader.o $(SRCDIR)/collisions.o $(SRCDIR)/statistics.o

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod *.mod
	rm -f $(EXEC)
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.mod
	@echo "Clean complete"

run: $(EXEC)
	cd $(RUNDIR) && ./monte_carlo_0d
