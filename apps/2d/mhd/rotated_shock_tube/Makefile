### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/mhd/lib/Makefile.defs

# overridden library objects
#
InitApp      = InitApp

ALL_OBJECTS = $(COMMON_OBJECTS) \
MHDParams.o \
InitialParams.o

ALL_SOURCES = $(COMMON_SOURCES) \
MHDParams.cpp \
InitialParams.cpp
 
TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/mhd/lib/ $(LIB_SYS_INCLUDE)


LOCAL_INCLUDE+=-I.
### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
