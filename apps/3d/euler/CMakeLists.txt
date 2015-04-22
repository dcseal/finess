append_current_app_name()

override(SOURCES AuxFunc.cpp DFluxFunc.cpp D2FluxFunc.cpp SourceTermFunc.cpp)

list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(cylindrical)
add_subdirectory(riemann_test)
add_subdirectory(smooth)
