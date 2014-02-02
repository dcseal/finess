
### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# compiler options -- see file $(FINESS)/config/Makefile.config

# overridden library objects
#
#include $(FINESS)/apps/semi_lagrangian/lib/2d/cart/Makefile.defs

## ------------------------------------------------------------------------- ##
## Choose which DogSolveUser method to try:
##
##   Options include: semi-lagrangian scheme, hybrid semi-lagrangian, and
##   two-derivative method.
##
#DogSolveUser        = $(FINESS)/lib/2d/cart/DogSolveTD
## ------------------------------------------------------------------------- ##


## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
#D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##


ALL_OBJECTS = \
$(COMMON_OBJECTS) \

ALL_SOURCES = \
$(COMMON_SOURCES) \

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
