### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/acoustics/lib/Makefile.defs

# compiler options
#
# overridden library objects
#
#FinSolveUser           = $(FINESS)/lib/1d/FinSolveRK_PIF
#ConstructIntegratedF    = ConstructIntegratedF
FinSolveUser            = $(APP_LIB)/FinSolveDiffTransform

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES)

### targets ###

include $(FINESS)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
