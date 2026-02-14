# Makefile for Monte Carlo Single Particle Trace Simulation
# Usage:
#   make trace    - Build single particle trace
#   make hist     - Build collision histogram
#   make all      - Build both
#   make run      - Run trace simulation
#   make run-hist - Run histogram analysis
#   make clean    - Remove build files

# Compiler settings
FC = gfortran
FFLAGS = -O2 -Wall -Wextra
FFLAGS_DEBUG = -O0 -g -fcheck=all -fbacktrace

# Directories
SRCDIR = code
RUNDIR = run

# Common source files
COMMON_SRCS = $(SRCDIR)/constants.f90 \
              $(SRCDIR)/types.f90 \
              $(SRCDIR)/cdf_reader.f90 \
              $(SRCDIR)/cross_sections.f90 \
              $(SRCDIR)/dynamics.f90

# Trace executable
TRACE_TARGET = $(RUNDIR)/trace_single_particle
TRACE_SRCS = $(COMMON_SRCS) $(SRCDIR)/trace_single_particle.f90

# Histogram executable
HIST_TARGET = $(RUNDIR)/collision_histogram
HIST_SRCS = $(SRCDIR)/constants.f90 \
            $(SRCDIR)/types.f90 \
            $(SRCDIR)/cdf_reader.f90 \
            $(SRCDIR)/cross_sections.f90 \
            $(SRCDIR)/collision_histogram.f90

# Module files
MODS = constants.mod types.mod cdf_reader.mod cross_sections.mod dynamics.mod

# Default target
all: trace hist

# Single particle trace
trace: $(TRACE_TARGET)

$(TRACE_TARGET): $(TRACE_SRCS)
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(TRACE_SRCS)
	@mv -f *.mod $(RUNDIR)/ 2>/dev/null || true
	@echo "Built: $@"

# Collision histogram
hist: $(HIST_TARGET)

$(HIST_TARGET): $(HIST_SRCS)
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(HIST_SRCS)
	@mv -f *.mod $(RUNDIR)/ 2>/dev/null || true
	@echo "Built: $@"

# Debug builds
debug-trace: FFLAGS = $(FFLAGS_DEBUG)
debug-trace: $(TRACE_TARGET)

debug-hist: FFLAGS = $(FFLAGS_DEBUG)
debug-hist: $(HIST_TARGET)

# Run simulations
run: $(TRACE_TARGET)
	cd $(RUNDIR) && ./trace_single_particle

run-hist: $(HIST_TARGET)
	cd $(RUNDIR) && ./collision_histogram

# Plot results
plot-trace: run
	cd $(RUNDIR) && gnuplot plot_trace.gp

plot-hist: run-hist
	cd $(RUNDIR) && gnuplot plot_histogram.gp

# Clean
clean:
	rm -f $(TRACE_TARGET) $(HIST_TARGET) *.mod $(RUNDIR)/*.mod

distclean: clean
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.png

.PHONY: all trace hist debug-trace debug-hist run run-hist plot-trace plot-hist clean distclean
