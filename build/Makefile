SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)
# The compiler
FC_PC=gfortran
# flags for debugging or for maximum performance, comment as necessary
LD_PC_FLAGS=-O3 -fopenmp -lgomp 
FC_PC_FLAGS=-O3 -g -c  -ffree-form -ffree-line-length-none -fopenmp -fbounds-check 

EXEC=SPEED2D

.PHONY: all
all: $(OBJS) $(EXEC)

$(EXEC): $(OBJS)
	$(FC_PC) -o $@  $(OBJS) $(LD_PC_FLAGS)

$(OBJS):%.o: %.f90 MODULES.o
	$(FC_PC) $(FC_PC_FLAGS) $< -o $@

.PHONY: clean
clean:
	-rm -f *.o *.mod SPEED2D
