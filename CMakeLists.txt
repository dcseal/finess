cmake_minimum_required (VERSION 2.8.3)
project(finess)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(override)

include(check_unix)
include(check_compiler)
include(check_python)
include(setup_silo)
include(setup_openmp)
include(concat_before_every_string_in_list)
include(copy_scripts)
include(add_generate_sources)
include(add_copy_ini)
include(append_current_app_name)

check_unix()
check_compiler()
check_python()

set(EXTERNAL_LIBS)

option(USE_SILO "Use Silo output format." ON)
if(USE_SILO)
    setup_silo()
endif()

option(USE_OPENMP "Use OpenMP." ON)
if(USE_OPENMP)
    setup_openmp()
endif()

copy_scripts()

configure_file(
    ${PROJECT_SOURCE_DIR}/lib/config/config.h.in
    ${PROJECT_BINARY_DIR}/lib/config/config.h
    )
include_directories(${PROJECT_BINARY_DIR}/lib/config)

include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/lib/sys)

set(SOURCES
	${PROJECT_SOURCE_DIR}/lib/dog_io.cpp
	${PROJECT_SOURCE_DIR}/lib/dog_math.cpp
	${PROJECT_SOURCE_DIR}/lib/dog_str.cpp
	${PROJECT_SOURCE_DIR}/lib/main_global.cpp
	${PROJECT_SOURCE_DIR}/lib/RunStartScript.cpp
	${PROJECT_SOURCE_DIR}/lib/sys/assert.cpp
	${PROJECT_SOURCE_DIR}/lib/sys/debug.cpp
	${PROJECT_SOURCE_DIR}/lib/tensors.cpp
	${PROJECT_SOURCE_DIR}/lib/WenoReconstruct.cpp
	)

set(TARGET_NAME_PREFIX)


add_copy_ini()

add_subdirectory(apps/1d)
add_subdirectory(apps/2d)
add_subdirectory(apps/3d)