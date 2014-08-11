### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/buckley_lev/lib/Makefile.defs

# overridden library objects
#
#DogSolveUser  = $(DOGPACK)/lib/1d/DogSolveRK4
ConSoln                     = ConSoln
ConstructIntegratedF        = ConstructIntegratedF

# place to add object modules
#
#ALL_OBJECTS = $(COMMON_OBJECTS) ...
#ALL_SOURCES = $(COMMON_SOURCES) ...

### targets ###

include $(DOGPACK)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
