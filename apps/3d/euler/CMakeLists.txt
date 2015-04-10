append_current_app_name()

list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/AuxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/D2FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/DFluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SourceTermFunc.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(cylindrical)
add_subdirectory(riemann_test)
add_subdirectory(smooth)
