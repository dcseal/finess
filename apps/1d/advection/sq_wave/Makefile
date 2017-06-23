### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/advection/lib/Makefile.defs

# compiler options
#
# overridden library objects
#
#ConstructIntegratedF    = ConstructIntegratedF
#FinSolveUser            = FinSolveUser
FinSolveUser            = $(APP_LIB)/FinSolveDiffTransform
#FinSolveUser            = $(APP_LIB)/FinSolveDiffTransformFluxLimiter

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES)

### targets ###

include $(FINESS)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
