# Makefile for Monte Carlo Neutral Particle Simulation
# Usage:
#   make        - Build the executable in run/
#   make clean  - Remove build files
#   make run    - Build and run

# Compiler settings
FC = gfortran
FFLAGS = -O2 -Wall -Wextra
FFLAGS_DEBUG = -O0 -g -fcheck=all -fbacktrace

# Directories
SRCDIR = code
RUNDIR = run

# Target executable (in run directory)
TARGET = $(RUNDIR)/monte_carlo

# Source files (order matters for dependencies)
SRCS = $(SRCDIR)/constants.f90 \
       $(SRCDIR)/types.f90 \
       $(SRCDIR)/cdf_reader.f90 \
       $(SRCDIR)/cross_sections.f90 \
       $(SRCDIR)/scoring.f90 \
       $(SRCDIR)/dynamics.f90 \
       $(SRCDIR)/io.f90 \
       $(SRCDIR)/main.f90

# Module files generated
MODS = constants.mod types.mod cdf_reader.mod cross_sections.mod \
       scoring.mod dynamics.mod io.mod

# Default target
all: $(TARGET)

# Build executable directly in run/
$(TARGET): $(SRCS)
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(SRCS)
	@mv -f *.mod $(RUNDIR)/ 2>/dev/null || true
	@echo "Built: $@"

# Debug build
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: $(TARGET)

# Run simulation
run: $(TARGET)
	cd $(RUNDIR) && ./monte_carlo

# Clean build files
clean:
	rm -f $(TARGET) $(MODS) *.mod $(RUNDIR)/*.mod

# Clean everything including results
distclean: clean
	rm -f $(RUNDIR)/results.csv

.PHONY: all debug run clean distclean
