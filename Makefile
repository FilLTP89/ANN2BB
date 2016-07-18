# SOURCE/BILD DIRECTORIES
BUILDDIR = build
SRCDIR   = speed2d_nl_beta/src
MODDIR   = speed2d_nl_beta/modules

# MACROS
# ---> Sources
SRCS=$(wildcard $(SRCDIR)/*.f90)
# ---> Modules
MODS=$(wildcard $(MODDIR)/*.f90)
# ---> Objects
OBJS   = $(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))
COBJS += $(BUILDDIR)/MODULES.o $(BUILDDIR)/SPEED_FIELDS.o
EOBJS += $(BUILDDIR)/WRITE_MONITORS.o $(BUILDDIR)/SEISMIC_FORCES.o $(BUILDDIR)/NONLINEAR.o
# ---> Executable
EXEC=$(BUILDDIR)/SPEED2D
# ---> Remove
RM:= rm -fr

# COMPILER 
FC_PC = gfortran
# COMPILATION FLAGS
LD_PC_FLAGS = -O3 -fopenmp -lgomp
FC_PC_FLAGS = -O3 -g -c -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -J build 

# DEPENDENCIES
$(EOBJS): $(COBJS)
$(OBJS) : $(EOBJS) $(COBJS)

# Make all instructions
.PHONY: all
all: $(COBJS) $(EOBJS) $(OBJS) $(EXEC) 

$(EXEC) : $(OBJS) $(COBJS) $(EOBJS)
	$(FC_PC) -o $@  $(OBJS) $(COBJS) $(EOBJS) $(LD_PC_FLAGS)

$(BUILDDIR)/%.o : $(SRCDIR)/%.f90 $(MODDIR)/%.f90 $(BUILDDIR)/WRITE_MONITORS.o $(BUILDDIR)/SPEED_FIELDS.f90 $(BUILDDIR)/NONLINEAR.o $(BUILDDIR)/SEISMIC_FORCES.o $(BUILDDIR)/MODULES.o
	$(FC_PC) $(FC_PC_FLAGS) $^ -o $@ 


# Make clean instructions
.PHONY: clean
clean:
	$(RM) $(BUILDDIR)/*.o $(BUILDDIR)/*.mod
