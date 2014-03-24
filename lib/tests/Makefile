include $(FINESS)/lib/Makefile.defs

DOGDEFS_SOURCES = \
 $(TOPLIB)/constants.h \
 $(LIB_SYS)/assert.h \
 $(LIB_SYS)/assert.cpp \
 $(LIB_SYS)/fcmp.c \
 $(LIB_SYS)/debug.h \
 $(LIB_SYS)/debug.cpp \
 $(TOPLIB)/tensors.h \
 $(TOPLIB)/tensors.cpp \

DOGDEFS_OBJECTS = \
 $(LIB_SYS)/assert.o \
 $(LIB_SYS)/fcmp.o \
 $(LIB_SYS)/debug.o \
 $(TOPLIB)/tensors.o \

DOGDEFS_DEPENDENCIES = \
 $(DOGDEFS_SOURCES) \
 $(DOGDEFS_OBJECTS) \

TEST_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test.o

TEST_CONSTANTS_OBJECTS = \
  $(CONSTANTS_OBJECTS) \
  test_constants.o

TEST_DOG_MATH_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test_dog_math.o

TEST_TENSOR_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  test_tensors.o

TIME_SOURCES = \
 $(DOGDEFS_SOURCES) \
 $(LIB_SYS)/ext_time.h \
 $(LIB_SYS)/ext_time.cpp \
 time_test.cpp \

TIME_OBJECTS = \
  $(DOGDEFS_OBJECTS) \
  $(LIB_SYS)/ext_time.o \
 #time_test.o \

TEST_WENORECONSTRUCT_OBJECTS = \
  test_wenoreconstruct.o

TEST_WENORECONSTRUCT_DEPENDENCIES = \
  $(TOPLIB)/tensors.o \
  $(TOPLIB)/WenoParams.o \
  $(TOPLIB)/DogParams.o  \
  $(TOPLIB)/sys/assert.o \
  $(TOPLIB)/sys/IniDocument.o \
  $(TOPLIB)/sys/debug.o  \
  $(TOPLIB)/dog_str.o  \
  $(TOPLIB)/dog_io.o  \
  $(TOPLIB)/dog_ini.o \
  $(TOPLIB)/WenoReconstruct.o



%.o : %.cpp
	$(LINK) $(FFLAGS) $< -o $@

#.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

EXECUTABLES = \
  test \
  test_tensors \
  test_constants \
  test_dog_math \
  time_test \
  test_wenoreconstruct

all: $(EXECUTABLES)

test: $(DOGDEFS_DEPENDENCIES) test.cpp test.o
	$(LINK) $(LFLAGS) $(TEST_OBJECTS) -o test

test_constants: $(DOGDEFS_DEPENDENCIES) test_constants.cpp test_constants.o
	$(LINK) $(LFLAGS) $(TEST_CONSTANTS_OBJECTS) -o test_constants

test_dog_math: $(DOGDEFS_DEPENDENCIES) test_dog_math.cpp test_dog_math.o
	$(LINK) $(LFLAGS) $(TEST_DOG_MATH_OBJECTS) -o test_dog_math

test_tensors:  $(DOGDEFS_DEPENDENCIES) test_tensors.cpp test_tensors.o
	$(LINK) $(LFLAGS) $(TEST_TENSOR_OBJECTS) -o test_tensors

time_test: time_test.o $(TIME_OBJECTS) $(TIME_SOURCES)
	$(COMPILE) $(TIME_FLAGS) time_test.cpp -o time_test.o
	$(LINK) time_test.o $(TIME_OBJECTS) -o time_test

test_wenoreconstruct: test_wenoreconstruct.o $(TEST_WENORECONSTRUCT_OBJECTS) $(TEST_WENORECONSTRUCT_DEPENDENCIES)
	$(LINK) $(LFLAGS) $(TEST_WENORECONSTRUCT_OBJECTS) $(TEST_WENORECONSTRUCT_DEPENDENCIES) -o test_wenoreconstruct

CLEAN_FILES = *.o $(EXECUTABLES)

include $(FINESS)/lib/Makefile.targets
