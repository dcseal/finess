# modifies variable: generated_sources
macro(add_generate_sources)
    set(ags_opt_main_python_script MAIN_PYTHON_SCRIPT)
    set(ags_opt_other_python_scripts OTHER_PYTHON_SCRIPTS)
    set(ags_opt_generated_sources GENERATED_SOURCES)
    cmake_parse_arguments(ags
                          ""
                          "${ags_opt_main_python_script}"
                          "${ags_opt_other_python_scripts};${ags_opt_generated_sources}"
                          ${ARGN})
    concat_before_every_string_in_list("${CMAKE_CURRENT_BINARY_DIR}/"
                                       ags_output_list
                                       ${ags_MAIN_PYTHON_SCRIPT}
                                       ${ags_OTHER_PYTHON_SCRIPTS}
                                       ${ags_GENERATED_SOURCES} )

    concat_before_every_string_in_list("${CMAKE_CURRENT_SOURCE_DIR}/"
                                       ags_scripts_list
                                       ${ags_MAIN_PYTHON_SCRIPT}
                                       ${ags_OTHER_PYTHON_SCRIPTS} )

    add_custom_command(OUTPUT ${ags_output_list}
                       # Very Unix-y.  Have yet to find more
                       # platform-indendendent solutions
                       COMMAND cp ${ags_scripts_list}
                                  ${CMAKE_CURRENT_BINARY_DIR}
                       COMMAND ${PROJECT_BINARY_DIR}/scripts/pythonpath_wrapper
                               ${PROJECT_BINARY_DIR}/python
                               ${PYTHON_EXECUTABLE}
                               ${ags_MAIN_PYTHON_SCRIPT}
                       DEPENDS ${ags_MAIN_PYTHON_SCRIPT}
                               ${ags_OTHER_PYTHON_SCRIPTS}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    
    concat_before_every_string_in_list(${CMAKE_CURRENT_BINARY_DIR}/
                                       generated_sources
                                       ${ags_GENERATED_SOURCES}  )
endmacro()
