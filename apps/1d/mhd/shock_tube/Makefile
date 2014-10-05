### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/1d/mhd/lib/Makefile.defs

# overridden library objects
#

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) 
ALL_SOURCES = $(COMMON_SOURCES)
TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/1d/mhd/lib/ $(LIB_SYS_INCLUDE)

### targets ###

include $(FINESS)/lib/1d/Makefile.targets

### DO NOT remove this line - make depends on it ###
