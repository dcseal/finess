# Add target copy_ini
# 	  which copies all the ini files from source tree to binary tree
# This target is not included in target 'all' by default

macro(add_copy_ini)
	add_custom_target(copy_ini
					  COMMAND ${PROJECT_SOURCE_DIR}/scripts/copy_ini
					  		  ${PROJECT_BINARY_DIR}
					  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endmacro()