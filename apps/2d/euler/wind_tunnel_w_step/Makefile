### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/euler/lib/Makefile.defs

# overridden library objects
#
InitApp      = InitApp

# Used to replace values inside wedge with (constant) initial conditions
AfterFullTimeStep   = AfterFullTimeStep

# Routines needed to be overwritten in order to implement the correct boundary
# conditions.  For this problem, we need two SetBndValues functions: one for
# x-direction, and one for y-direction
ConstructL          = ConstructL

ALL_OBJECTS = $(COMMON_OBJECTS) \
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
EulerParams.cpp \

TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/euler/lib/ $(LIB_SYS_INCLUDE)

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
