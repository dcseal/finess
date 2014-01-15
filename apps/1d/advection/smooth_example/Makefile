### definitions ###

# default/library definitions
include $(FINESS)/lib/1d/Makefile.defs

# compiler options
#
# overridden library objects
#

## ------------------------------------------------------------------------- ##
# List of a few different options explored for this problem:
#
#   DogSolveRK4  : classical, fourth-order time integration.
#   DogSolveSL   : Semi-Lagrangian solver
#   DogSolveTDRK : Fourth-order, two-derivative solver
#
# We currently have no mechanism in place to select out of these various
# solvers, so this Makefile has to link to the appropriate solver.
#
#DogSolveUser           = $(DOGPACK)/apps/semi_lagrangian/lib/1d/DogSolveSL
#DogSolveUser           = $(DOGPACK)/lib/1d/DogSolveRK4
#DogSolveUser           = $(DOGPACK)/lib/1d/DogSolveTDRK
## ------------------------------------------------------------------------- ##

#RiemannSolveLxW         = $(DOGPACK)/lib/1d/LaxWendroff/RiemannSolveLxW_HLL

## ------------------------------------------------------------------------- ##
# These time derivatives are required for the Lax-Wendroff (and
# two-derivative) code:
DFluxFunc              = DFluxFunc
D2FluxFunc             = D2FluxFunc
## ------------------------------------------------------------------------- ##

# place to add object modules
#
ALL_OBJECTS = $(COMMON_OBJECTS) $(SL_LIB_OBJECTS)
ALL_SOURCES = $(COMMON_SOURCES) $(SL_LIB_SOURCES)

### targets ###

include $(DOGPACK)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
