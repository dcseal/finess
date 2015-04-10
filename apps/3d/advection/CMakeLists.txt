append_current_app_name()

list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/D2FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/DFluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(smooth_example)
