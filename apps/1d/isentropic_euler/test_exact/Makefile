### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/isentropic_euler/lib/Makefile.defs

# overridden library objects
#
FinSolveUser            = $(APP_LIB)/FinSolveDiffTransform


# uncomment these lines to use the HLL(E) Riemann solver
#RiemannSolve     = $(LIB1D)/RiemannSolveHLL
#RiemannSolveLxW  = $(LIB1D)/LaxWendroff/RiemannSolveLxW_HLL

# place to add object modules
#
#ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS)
#ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES)

ALL_OBJECTS = $(COMMON_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES)

TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/1d/isentropic_euler/lib/ $(LIB_SYS_INCLUDE)

### targets ###

include $(FINESS)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
