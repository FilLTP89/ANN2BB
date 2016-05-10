# Source/Build Folders and targets
SRCDIR=src
BUILDDIR=build
# Macros
SRCS=$(wildcard $(SRCDIR)/*.f90)
OBJS=$(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))
EXEC=SPEED2D
# The compiler
FC_PC=gfortran
# flags for debugging or for maximum performance, comment as necessary
LD_PC_FLAGS= -fopenmp -lgomp
FC_PC_FLAGS= -g -c  -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -J build 
# Make all instructions
.PHONY: all
all: dir $(BUILDDIR)/$(EXEC) 

dir:
	-mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXEC): $(OBJS)
	$(FC_PC) -o $@ $(OBJS) $(LD_PC_FLAGS) 

$(OBJS) : $(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(BUILDDIR)/MODULES.o
	$(FC_PC) $(FC_PC_FLAGS) $< -o $@

# Make clean instructions
.PHONY: clean
clean:
	-rm -fr $(BUILDDIR) $(SRCDIR)/*.mod
