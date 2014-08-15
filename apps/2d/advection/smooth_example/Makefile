### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/advection/lib/Makefile.defs

# compiler options -- see file $(FINESS)/config/Makefile.config

## ------------------------------------------------------------------------- ##
# overridden library objects
#
AuxFunc                 = AuxFunc
DogSolveUser            = $(FINESS)/lib/2d/FinSolveRK_PIF
## ------------------------------------------------------------------------- ##

ALL_OBJECTS = \
$(COMMON_OBJECTS) \

ALL_SOURCES = \
$(COMMON_SOURCES) \

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
