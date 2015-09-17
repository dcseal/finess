### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/euler/lib/Makefile.defs

# overridden library objects
#

# Used to replace values inside wedge with (constant) initial conditions
AfterFullTimeStep   = AfterFullTimeStep

# Routines needed to be overwritten in order to implement the correct boundary
# conditions.  For this problem, we need two SetBndValues functions: one for
# x-direction, and one for y-direction
ConstructL          = ConstructL
ConstructLxWL	    = ConstructLxWL
FinSolveLxW	    	= FinSolveLxW

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS) \

ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES) \

TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/euler/lib/ $(LIB_SYS_INCLUDE)

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
