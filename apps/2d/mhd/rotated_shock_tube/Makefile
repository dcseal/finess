### definitions ###

# default/library definitions
include $(FINESS)/lib/2d/Makefile.defs

# Default (application level) definintions
include $(FINESS)/apps/2d/mhd/lib/Makefile.defs

# overridden library objects
#
#InitApp      = InitApp

ALL_OBJECTS = $(COMMON_OBJECTS) \
MHDParams.o \
InitialParams.o \
IniParams.o

ALL_SOURCES = $(COMMON_SOURCES) \
MHDParams.cpp \
InitialParams.cpp \
IniParams.cpp


TOPLIB_INCLUDE = -I$(TOPLIB) -I$(FINESS)/apps/2d/mhd/lib/ $(LIB_SYS_INCLUDE)


LOCAL_INCLUDE+=-I.




### targets ###

include $(FINESS)/lib/2d/Makefile.targets

IniParams.cpp: generate_iniparams.py
	python $^

IniParser.h: IniParser.rl
	ragel $^ -o $@
 
### DO NOT remove this line - make depends on it ###
