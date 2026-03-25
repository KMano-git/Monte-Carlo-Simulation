# Makefile for 1D3V Monte Carlo Simulation
#===============================================================================

FC = gfortran
FFLAGS = -O2 -Wall -Wextra -fcheck=bounds -J$(MODDIR)
LDFLAGS =

ifdef OPENMP
  FFLAGS += -fopenmp
  LDFLAGS += -fopenmp
endif

SRCDIR = code
RUNDIR = run
MODDIR = module

SOURCES = $(SRCDIR)/constants.f90 \
          $(SRCDIR)/data_types.f90 \
          $(SRCDIR)/random_utils.f90 \
          $(SRCDIR)/cdf_reader.f90 \
          $(SRCDIR)/cross_sections.f90 \
          $(SRCDIR)/dynamics.f90 \
          $(SRCDIR)/scoring.f90 \
          $(SRCDIR)/io.f90 \
          $(SRCDIR)/main.f90

OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(MODDIR)/%.o,$(SOURCES))
EXEC = $(RUNDIR)/monte_carlo_1d3v

.PHONY: all clean run

all: $(EXEC)

$(EXEC): $(OBJECTS)
	@mkdir -p $(RUNDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Build complete: $@"

$(MODDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(MODDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(MODDIR)/data_types.o: $(MODDIR)/constants.o
$(MODDIR)/random_utils.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o
$(MODDIR)/cdf_reader.o: $(MODDIR)/constants.o
$(MODDIR)/cross_sections.o: $(MODDIR)/constants.o $(MODDIR)/cdf_reader.o
$(MODDIR)/dynamics.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
	$(MODDIR)/random_utils.o $(MODDIR)/cross_sections.o $(MODDIR)/cdf_reader.o
$(MODDIR)/scoring.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o \
	$(MODDIR)/random_utils.o $(MODDIR)/cross_sections.o $(MODDIR)/cdf_reader.o
$(MODDIR)/io.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o $(MODDIR)/random_utils.o
$(MODDIR)/main.o: $(MODDIR)/constants.o $(MODDIR)/data_types.o $(MODDIR)/random_utils.o \
	$(MODDIR)/cdf_reader.o $(MODDIR)/cross_sections.o $(MODDIR)/dynamics.o \
	$(MODDIR)/scoring.o $(MODDIR)/io.o

clean:
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod
	rm -f $(EXEC)
	rm -f $(RUNDIR)/monte_carlo
	rm -f $(RUNDIR)/*.csv $(RUNDIR)/*.dat $(RUNDIR)/*.mod $(RUNDIR)/*.png
	rm -rf $(RUNDIR)/figure
	@echo "Clean complete"

run: $(EXEC)
	cd $(RUNDIR) && ./monte_carlo_1d3v
