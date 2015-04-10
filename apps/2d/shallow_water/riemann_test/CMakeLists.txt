append_current_app_name()

list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/QinitFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetBndValues.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SourceTermFunc.cpp
	)

add_generate_sources(MAIN_PYTHON_SCRIPT generate_iniparams.py
					GENERATED_SOURCES IniParams.h IniParams.cpp)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(TARGET_NAME ${TARGET_NAME_PREFIX}_finess.exe)
add_executable(${TARGET_NAME} ${SOURCES} ${generated_sources})
set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME finess.exe)
target_link_libraries(${TARGET_NAME} ${EXTERNAL_LIBS})