append_current_app_name()

override(SOURCES AuxFunc.cpp DFluxFunc.cpp D2FluxFunc.cpp)
override(SOURCES ApplyMPPLimiter2D.cpp)

list(APPEND SOURCES
#	${CMAKE_CURRENT_SOURCE_DIR}/ApplyMPPLimiter2D.cpp
#	${CMAKE_CURRENT_SOURCE_DIR}/AuxFunc.cpp
#	${CMAKE_CURRENT_SOURCE_DIR}/D2FluxFunc.cpp
#	${CMAKE_CURRENT_SOURCE_DIR}/DFluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(mach_reflection)
add_subdirectory(radial_shock)
add_subdirectory(riemann_test)
add_subdirectory(sedov_blast)
add_subdirectory(shock_diffraction)
add_subdirectory(test_exact1)
add_subdirectory(test_exact2)
add_subdirectory(wind_tunnel_w_step)