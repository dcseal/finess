### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# compiler options -- see file $(FINESS)/config/Makefile.config

## ------------------------------------------------------------------------- ##
# overridden library objects
#

# These time derivatives are required for the Lax-Wendroff code:
DFluxFunc              = DFluxFunc
D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##

ALL_OBJECTS = \
$(COMMON_OBJECTS) \

ALL_SOURCES = \
$(COMMON_SOURCES) \

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
