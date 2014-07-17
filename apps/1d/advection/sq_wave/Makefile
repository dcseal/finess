### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/advection/lib/Makefile.defs

# compiler options
#
# overridden library objects
#
ConstructIntegratedF    = ConstructIntegratedF

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(SL_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(SL_LIB_SOURCES)

### targets ###

include $(FINESS)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
