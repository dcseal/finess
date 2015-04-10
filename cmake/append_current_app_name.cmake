macro(append_current_app_name)
	get_filename_component(current_app_name ${CMAKE_CURRENT_SOURCE_DIR} NAME)
	set(TARGET_NAME_PREFIX "${TARGET_NAME_PREFIX}_${current_app_name}")
endmacro()