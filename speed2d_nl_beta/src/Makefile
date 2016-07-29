SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

all: $(OBJS)

clean:
	$(RM) $(OBJS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
