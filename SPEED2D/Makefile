# SOURCE/BILD DIRECTORIES
BUILDDIR = speed2d_nl_beta/build
SRCDIR   = speed2d_nl_beta/src

# MACROS
# ---> Sources
SRCS=$(wildcard $(SRCDIR)/*.f90)
# ---> Objects
OBJS=$(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(SRCS))
# ---> Executable
EXEC=$(BUILDDIR)/SPEED2D
# ---> Remove
RM:= rm -fr


# COMPILER 
FC_PC=gfortran
# COMPILATION FLAGS
LD_PC_FLAGS=-O3 -fopenmp -lgomp
FC_PC_FLAGS=-O3 -g -c -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -J speed2d_nl_beta/build 

# Make all instructions
.PHONY: all
all: $(OBJS) $(EXEC)


$(EXEC): $(OBJS) 
	$(FC_PC) -o $@ $(OBJS) $(LD_PC_FLAGS)
	
$(OBJS): $(BUILDDIR)/%.o: $(SRCDIR)/%.f90 $(BUILDDIR)/MODULES.o $(BUILDDIR)/WRITE_MONITORS.o $(BUILDDIR)/SEISMIC_FORCES.o $(BUILDDIR)/SPEED_FIELDS.o
	$(FC_PC) $(FC_PC_FLAGS)  $< -o $@

.PHONY: srcdir
srcdir:
	@echo $(SRCDIR)

.PHONY: src
src:
	@echo $(SRCS)

.PHONY: obj
obj:
	@echo $(OBJS)

# Make clean instructions
.PHONY: clean
clean:
	$(RM) $(BUILDDIR)/*.o $(BUILDDIR)/*.mod
