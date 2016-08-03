SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

all: $(OBJS)

clean:
	$(RM) $(OBJS) *.mod

SEISMIC_FORCES.o : SEISMIC_FORCES.f90 SPEED_FIELDS.o
	$(FC) $(FCFLAGS) -c -o $@ $<

NONLINEAR.o : NONLINEAR.f90 SPEED_FIELDS.o
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
