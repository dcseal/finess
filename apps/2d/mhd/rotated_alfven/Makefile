### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/mhd/lib/Makefile.defs

# overridden library objects
FinSolveRK = FinSolveRK
FinSolveLxW = FinSolveLxW

AfterQinit = AfterQinit
UpdateSoln = UpdateSoln
AfterStep = AfterStep
GetCFL = GetCFL

# Hamilton-Jacobi rhs for magnetic potential
ConstructHJ_L = ConstructHJ_L

ALL_OBJECTS = $(COMMON_OBJECTS) \
$(ConstructHJ_L).o

ALL_SOURCES = $(COMMON_SOURCES) \
$(ConstructHJ_L).cpp \
$(ConstructHJ_L).h



TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/mhd/lib/ $(LIB_SYS_INCLUDE)


### targets ###

include $(FINESS)/lib/2d/Makefile.targets



### DO NOT remove this line - make depends on it ###
