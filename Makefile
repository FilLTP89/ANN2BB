# Source/Build Folders and targets
SRCDIR=speed2d_nl_beta/src
BUILDDIR=build
MODDIR=speed2d_nl_beta/modules

# Macros
SRCS=$(wildcard $(SRCDIR)/*.f90)
MODS=$(wildcard $(MODDIR)/*.f90)
OBJS=$(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))
COBJS += $(MODDIR)/MODULES.o $(MODDIR)/SPEED_FIELDS.o
EOBJS += $(MODDIR)/WRITE_MONITORS.o $(MODDIR)/SEISMIC_FORCES.o $(MODDIR)/NONLINEAR.o
# Executable
EXEC=$(BUILDDIR)/SPEED2D

# Remove
RM:= rm -fr

# The compiler
FC_PC = gfortran
LD_PC_FLAGS = -O3 -fopenmp -lgomp
FC_PC_FLAGS = -O3 -g -c -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -J build 

# Dependencies
$(EOBJS) : $(COBJS)
$(OBJS)  : $(EOBJS) $(COBJS) 

# Make all instructions
.PHONY: all
all: dir $(OBJS) $(EXEC) 

dir:
	-mkdir -p $(BUILDDIR)

$(BUILDDIR)/SPEED_FIELDS.o : $(MODDIR)/SPEED_FIELD.f90
	@echo '$@'	
	$(FC_PC) $(FC_PC_FLAGS) -o $<
	@echo 'OBJECT' 
	@echo '$@'

$(BUILDDIR)/MODULES.o : $(MODDIR)/MODULES.f90
	@echo '$@'	
	$(FC_PC) $(FC_PC_FLAGS) -o $<
	@echo 'OBJECT' 
	@echo '$@'

$(BUILDDIR)/NONLINEAR.o : $(MODDIR)/NONLINEAR.f90
	@echo '$@'	
	$(FC_PC) $(FC_PC_FLAGS) -o $<
	@echo 'OBJECT' 
	@echo '$@'

$(BUILDDIR)/WRITE_MONITORS.o : $(MODDIR)/WRITE_MONITORS.f90
	@echo '$@'	
	$(FC_PC) $(FC_PC_FLAGS) -o $<
	@echo 'OBJECT' 
	@echo '$@'

$(BUILDDIR)/SEISMIC_FORCES.o : $(MODDIR)/SEISMIC_FORCES.f90
	@echo '$@'	
	$(FC_PC) $(FC_PC_FLAGS) -o $<
	@echo 'OBJECT' 
	@echo '$@'

$(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(EOBJS) $(COBJS) 
	$(FC_PC) $(FC_PC_FLAGS) -o $^
	@echo 'OBJECT' 
	@echo '$@'

$(EXEC) : $(OBJS)
	@echo 'EXECUTABLE' 
	@echo '$@'
	$(FC_PC) -o $@ $^ $(LD_PC_FLAGS)
	@echo 'TARGET'
	@echo '$<'

# Make clean instructions
.PHONY: clean
clean:
	$(RM) $(BUILDDIR)/*.o $(BUILDDIR)/*.mod
