# copy the following directories to binary directory:
#	python
#	scripts
#	util
#	viz

macro(copy_scripts)
	add_custom_target(copy_python
						ALL	${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/python ${PROJECT_BINARY_DIR}/python
					   )
	add_custom_target(copy_scripts
						ALL ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/scripts ${PROJECT_BINARY_DIR}/scripts
					   )
	add_custom_target(copy_util
						ALL ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/util ${PROJECT_BINARY_DIR}/util
					   )
	add_custom_target(copy_viz
						ALL ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/viz ${PROJECT_BINARY_DIR}/viz
					   )

endmacro()