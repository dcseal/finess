append_current_app_name()

override(SOURCES ApplyMPPLimiter2D.cpp AfterStep.cpp FinSolveLxW.cpp FinSolveRK.cpp GetCFL.cpp UpdateSoln.cpp)
override(SOURCES DFluxFunc.cpp D2FluxFunc.cpp)
override(SOURCES AuxFunc.cpp)

list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/ApplyMPPLimiter2D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ConstructHJ_L.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(blast)
add_subdirectory(cloud_shock)
add_subdirectory(orszag_tang)
add_subdirectory(rotated_alfven)
add_subdirectory(rotated_shock_tube)
add_subdirectory(shock_tube)
add_subdirectory(vacuum_vortex)
add_subdirectory(vortex)