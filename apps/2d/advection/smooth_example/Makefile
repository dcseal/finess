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
FinSolveUser            = $(APP_LIB)/FinSolveDiffTransform
## ------------------------------------------------------------------------- ##

ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS)

ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES)


LOCAL_INCLUDE+=-I.
LOCAL_INCLUDE+=-I$(FINESS)/apps/2d/advection/lib/

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
