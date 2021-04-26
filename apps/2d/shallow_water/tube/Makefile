### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/shallow_water/lib/Makefile.defs

## ------------------------------------------------------------------------- ##
# overridden library objects
#
#Output                 = Output
DogSolveUser            = $(FINESS)/lib/2d/FinSolveRK_PIF
## ------------------------------------------------------------------------- ##

# place to add object modules
#
ALL_OBJECTS = \
$(COMMON_OBJECTS) \

ALL_SOURCES = \
$(COMMON_SOURCES) \

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
