append_current_app_name()

override(SOURCES FinSolveLxW.cpp AfterStep.cpp 
		 AuxFunc.cpp D2FluxFunc.cpp DFluxFunc.cpp)

list(APPEND SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/ApplyMPPLimiter3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ConstructHJ_L.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FluxFunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectLeftEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ProjectRightEig.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SetWaveSpd.cpp
	)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(blast)
add_subdirectory(cloud_shock)
add_subdirectory(rotated_alfven)
add_subdirectory(shocktubeX)
add_subdirectory(shocktubeY)
add_subdirectory(shocktubeZ)
add_subdirectory(smoothX)
add_subdirectory(smoothY)
add_subdirectory(smoothZ)
