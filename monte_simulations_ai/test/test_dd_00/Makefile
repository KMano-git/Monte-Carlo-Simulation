FC = gfortran
FFLAGS = -O3 -Wall

SRCDIR = ../../3d3v_event-driven/code

OBJS = constants.o cdf_reader.o main.o

all: main

constants.o: $(SRCDIR)/constants.f90
	$(FC) $(FFLAGS) -c $<

cdf_reader.o: $(SRCDIR)/cdf_reader.f90 constants.o
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 cdf_reader.o constants.o
	$(FC) $(FFLAGS) -c $<

main: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

clean:
	rm -f *.o *.mod main dcs_output.csv dcs_plot.png
