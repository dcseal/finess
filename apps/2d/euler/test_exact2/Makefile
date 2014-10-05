### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/euler/lib/Makefile.defs

# overridden library objects
#
SetWaveSpd      = SetWaveSpd

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(APP_LIB_OBJECTS) \
#EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) $(APP_LIB_SOURCES) \
#EulerParams.cpp \

TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/euler/lib/ $(LIB_SYS_INCLUDE)

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
