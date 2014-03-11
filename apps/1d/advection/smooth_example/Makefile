### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# compiler options
#
# overridden library objects
#
# DogSolveUser           = DogSolveUser

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(SL_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(SL_LIB_SOURCES)

### targets ###

include $(FINESS)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
