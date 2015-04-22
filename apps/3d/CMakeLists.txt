set(TARGET_NAME_PREFIX "apps_3d")
list(APPEND SOURCES
	${PROJECT_SOURCE_DIR}/lib/3d/blanks/AuxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/blanks/SourceTermFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/LaxWendroff/DFluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/LaxWendroff/D2FluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/AfterFullTimeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/AfterQinit.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/AfterStep.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/BeforeFullTimeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/BeforeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/ConSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/ConstructIntegratedR.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/ConstructL.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/ConstructLFL.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/ConstructLxWL.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/FinSolveLxW.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/FinSolveRK.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/FinSolveUser.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/GetCFL.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/GlobalWaveSpd.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/Output.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/RunFinpack.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/SampleFunction.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/SetRKinfo.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/UpdateSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/3d/WriteQhelp.cpp
)

include_directories(${PROJECT_SOURCE_DIR}/lib/3d)

add_subdirectory(advection)
add_subdirectory(burgers)
add_subdirectory(mhd)
add_subdirectory(euler)
