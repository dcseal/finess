append_current_app_name()

override(SOURCES ApplyMPPLimiter1D.cpp)

override(SOURCES AuxFunc.cpp DFluxFunc.cpp D2FluxFunc.cpp)
list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(blast)
add_subdirectory(harten_shock_tube)
add_subdirectory(sedov_blast)
add_subdirectory(shock_entropy)
add_subdirectory(shock_tube)
add_subdirectory(test_exact)