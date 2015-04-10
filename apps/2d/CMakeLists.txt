set(TARGET_NAME_PREFIX "apps_2d")
list(APPEND SOURCES
	${PROJECT_SOURCE_DIR}/lib/2d/AfterFullTimeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/AfterQinit.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/AfterStep.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/AfterUpdateSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ApplyLimiter.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ApplyMPPLimiter2D.cpp
	#${PROJECT_SOURCE_DIR}/lib/2d/ApplyPosLimiter.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/BeforeFullTimeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/BeforeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/blanks/AuxFunc.cpp
	# ${PROJECT_SOURCE_DIR}/lib/2d/blanks/FluxFunc.cpp
	# ${PROJECT_SOURCE_DIR}/lib/2d/blanks/ProjectLeftEig.cpp
	# ${PROJECT_SOURCE_DIR}/lib/2d/blanks/ProjectRightEig.cpp
	# ${PROJECT_SOURCE_DIR}/lib/2d/blanks/QinitFunc.cpp
	# ${PROJECT_SOURCE_DIR}/lib/2d/blanks/SetBndValues.cpp
	# ${PROJECT_SOURCE_DIR}/lib/2d/blanks/SetWaveSpd.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/blanks/SourceTermFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConstructIntegratedR.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConstructL.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConstructLFL.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConstructLxWL.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConvertQtoW.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ConvertWtoQ.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/FinSolveLxW.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/FinSolveMD.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/FinSolveRK.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/FinSolveRK_PIF.cpp
	#${PROJECT_SOURCE_DIR}/lib/2d/FinSolveUser.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/GetCFL.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/GlobalWaveSpd.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/GridSetup.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/LaxWendroff/D2FluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/LaxWendroff/DFluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/Output.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/ResInt.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/RunFinpack.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/SampleFunction.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/SetRKinfo.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/UpdateSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/2d/WriteQhelp.cpp
)

include_directories(${PROJECT_SOURCE_DIR}/lib/2d)

add_subdirectory(advection)
add_subdirectory(mhd)
add_subdirectory(euler)
add_subdirectory(shallow_water)