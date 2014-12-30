### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/mhd/lib/Makefile.defs

# overridden library objects
AfterQinit = AfterQinit
# Hamilton-Jacobi rhs for magnetic potential

ALL_OBJECTS = $(COMMON_OBJECTS) 

ALL_SOURCES = $(COMMON_SOURCES)

TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/mhd/lib/ $(LIB_SYS_INCLUDE)


### targets ###

include $(FINESS)/lib/2d/Makefile.targets



### DO NOT remove this line - make depends on it ###
