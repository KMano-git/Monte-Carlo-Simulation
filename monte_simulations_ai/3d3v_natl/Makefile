# Makefile for 3D3V Non-Analog Track-Length Monte Carlo Simulation
#===============================================================================

# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -fcheck=bounds -J$(MODDIR)
LDFLAGS = 

# Directories
SRCDIR = code
RUNDIR = run
MODDIR = module

# Source files (in compilation order - dependencies first)
SOURCES = $(SRCDIR)/constants.f90 \
          $(SRCDIR)/data_types.f90 \
          $(SRCDIR)/random_utils.f90 \
          $(SRCDIR)/cdf_reader.f90 \
          $(SRCDIR)/cross_sections.f90 \
          $(SRCDIR)/scoring.f90 \
          $(SRCDIR)/dynamics.f90 \
          $(SRCDIR)/io.f90 \
          $(SRCDIR)/main.f90

# Object files
OBJECTS = $(patsubst $(SRCDIR)/%.f90, $(MODDIR)/%.o, $(SOURCES))

# Executable
EXEC = $(RUNDIR)/monte_carlo_3d3v_natl

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
$(MODDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(MODDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies
$(MODDIR)/data_types.o: $(MODDIR)/constants.o
$(MODDIR)/random_utils.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o
$(MODDIR)/cdf_reader.o: $(MODDIR)/constants.o
$(MODDIR)/cross_sections.o: $(MODDIR)/constants.o $(MODDIR)/cdf_reader.o
$(MODDIR)/scoring.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
                     $(MODDIR)/cross_sections.o $(MODDIR)/random_utils.o $(MODDIR)/cdf_reader.o
$(MODDIR)/dynamics.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
                      $(MODDIR)/random_utils.o $(MODDIR)/cross_sections.o $(MODDIR)/cdf_reader.o
$(MODDIR)/io.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o $(MODDIR)/random_utils.o
$(MODDIR)/main.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o $(MODDIR)/random_utils.o \
                  $(MODDIR)/cdf_reader.o $(MODDIR)/cross_sections.o \
                  $(MODDIR)/scoring.o $(MODDIR)/dynamics.o $(MODDIR)/io.o

clean:
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod
	rm -f $(EXEC)
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.mod
	@echo "Clean complete"

run: $(EXEC)
	cd $(RUNDIR) && ./monte_carlo_3d3v_natl
