# Source/Build Folders and targets
SRCDIR=src
BUILDDIR=build
# Macros
SRCS=$(wildcard $(SRCDIR)/*.f90)
CSRCS += $(SRCDIR)/SPEED_FIELDS.f90 $(SRCDIR)/MODULES.f90
ESRCS += $(SRCDIR)/WRITE_MONITORS.90 $(SRCDIR)/NONLINEAR.f90 $(SRCDIR)/SEISMIC_FORCES.f90  
OBJS=$(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))
COBJS += $(BUILDDIR)/SPEED_FIELDS.o $(BUILDDIR)/MODULES.o
EOBJS += $(BUILDDIR)/WRITE_MONITORS.o $(BUILDDIR)/NONLINEAR.o $(BUILDDIR)/SEISMIC_FORCES.o
# Executable
EXEC=SPEED2D
# Remove
RM:= rm -fr
# The compiler
FC_PC = gfortran
LD_PC_FLAGS = -O3 -fopenmp -lgomp -J build
FC_PC_FLAGS = -O3 -g -c -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -J build 

# Dependencies
$(EOBJS) : $(COBJS)

$(OBJS) : $(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(FC_PC) -o $@ $(FC_PC_FLAGS) $^
	@echo 'OBJECT' 
	@echo '$@' 
	@echo 'TARGET'
	@echo '$^'

# Make all instructions
.PHONY: all
all: dir $(OBJS) $(BUILDDIR)/$(EXEC) 

dir:
	-mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(EXEC) : $(OBJS)
	@echo 'EXECUTABLE' 
	@echo '$@'
	$(FC_PC) -o $@  $(OBJS) $(LD_PC_FLAGS)
	@echo 'TARGET'
	@echo '$<'

# Make clean instructions
.PHONY: clean
clean:
	$(RM) $(BUILDDIR)/*.o $(BUILDDIR)/*.mod
