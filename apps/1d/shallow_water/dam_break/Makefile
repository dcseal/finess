
### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/shallow_water/lib/Makefile.defs

# overridden library objects
#
FinSolveUser           = $(APP_LIB)/FinSolveDiffTransform

# Experimental new method that uses WENO *differentiation* on the interior of
# the time-integrated flux function.
#ConstructIntegratedF    = ConstructIntegratedF

# place to add object modules
#
#ALL_OBJECTS = $(COMMON_OBJECTS) ...
#ALL_SOURCES = $(COMMON_SOURCES) ...

ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES)

### targets ###

include $(FINESS)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
