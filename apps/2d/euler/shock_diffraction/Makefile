### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

### overridden library objects ###

# For saving 'eulerhelp.dat', which contains gamma gas constant
AfterQinit          = AfterQinit

# Used to call EulerParams (saves gamma gas constant)
InitApp             = InitApp

# Lax-Wendroff routines
DFluxFunc           = DFluxFunc
D2FluxFunc          = D2FluxFunc

# Used to replace values inside wedge with (constant) initial conditions
AfterFullTimeStep   = AfterFullTimeStep

# Routines needed to be overwritten in order to implement the correct boundary
# conditions.  For this problem, we need two SetBndValues functions: one for
# x-direction, and one for y-direction
# ConstructL          = ConstructL

ALL_OBJECTS = $(COMMON_OBJECTS) \
EulerParams.o \

ALL_SOURCES = $(COMMON_SOURCES) \
EulerParams.cpp \

### targets ###

include $(FINESS)/lib/2d/Makefile.targets

### DO NOT remove this line - make depends on it ###
