set(TARGET_NAME_PREFIX "apps_1d")
list(APPEND SOURCES
	${PROJECT_SOURCE_DIR}/lib/1d/GetCFL.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/FinSolveRK.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConstructL.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/TimeStepSDC.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/UpdateSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/FinSolveMD.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ResInt.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/LaxWendroff/DFluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/LaxWendroff/D2FluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConstructIntegratedF.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/BeforeFullTimeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/AfterUpdateSoln.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/SetRKinfo.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/RunFinpack.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/BeforeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/AfterQinit.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/Output.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConstructLFL.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/FinSolveSDC.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/CopyQ.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/blanks/SetWaveSpd.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/blanks/AuxFunc.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/blanks/QinitFunc.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/blanks/ProjectLeftEig.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/blanks/ProjectRightEig.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/blanks/SetBndValues.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/blanks/SourceTermFunc.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/blanks/FluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/GlobalWaveSpd.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConvertQtoW.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/WriteQhelp.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/main.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/SetWaveSpd.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/AuxFunc.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/Output.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/QinitFunc.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/ProjectLeftEig.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/ProjectRightEig.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/SetBndValues.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/SourceTermFunc.cpp
	# ${PROJECT_SOURCE_DIR}/lib/1d/tests/weno_diff_order/FluxFunc.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/FinSolveLxW.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/StepSDC.cpp
	#${PROJECT_SOURCE_DIR}/lib/1d/FinSolveUser.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/Output_Extra.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/AfterFullTimeStep.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConvertWtoQ.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ConstructLxWL.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/FinSolveRK_PIF.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/ApplyMPPLimiter1D.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/SampleFunction.cpp
	${PROJECT_SOURCE_DIR}/lib/1d/AfterStep.cpp
)

include_directories(${PROJECT_SOURCE_DIR}/lib/1d)

add_subdirectory(advection)
add_subdirectory(burgers)
add_subdirectory(buckley_lev)
add_subdirectory(euler)
add_subdirectory(euler_radial_of_spherically_symmetric_3d)
add_subdirectory(mhd)

