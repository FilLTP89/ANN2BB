# Source/Build Folders and targets
SRCDIR=src
BUILDDIR=build
# Macros
SRCS=$(wildcard $(SRCDIR)/*.f90)
OBJS=$(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))
EXEC=SPEED2D
RM:= rm -rf
# The compiler
FC_PC=gfortran
# flags for debugging or for maximum performance, comment as necessary
LD_PC_FLAGS=-O3 -fopenmp -lgomp
FC_PC_FLAGS=-O3 -g -c  -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -J build 

# Dependencies
MODULES.o           : SPEED_FIELDS.o
NONLINEAR.o         : SPEED_FIELDS.o
WRITE_MONITORS.o    : SPEED_FIELDS.o

# Make all instructions
.PHONY: all
all: dir $(BUILDDIR)/$(EXEC) 

dir:
	-mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXEC): $(OBJS)
	$(FC_PC) -o $@ $(OBJS) $(LD_PC_FLAGS) 

$(OBJS) : $(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(BUILDDIR)/SPEED_FIELDS.o
	$(FC_PC) $(FC_PC_FLAGS) $< -o $@

$(OBJS) : $(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(BUILDDIR)/MODULES.o
	$(FC_PC) $(FC_PC_FLAGS) $< -o $@

$(OBJS) : $(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(BUILDDIR)/NONLINEAR.o
	$(FC_PC) $(FC_PC_FLAGS) $< -o $@

$(OBJS) : $(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(BUILDDIR)/WRITE_MONITORS.o
	$(FC_PC) $(FC_PC_FLAGS) $< -o $@

# Make clean instructions
.PHONY: clean
clean:
	-rm -fr $(BUILDDIR) $(SRCDIR)/*.mod
