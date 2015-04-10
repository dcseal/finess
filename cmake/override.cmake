function(override sources_list)
    foreach (original_file ${${sources_list}})
        set(overridden_file ${original_file})
        foreach(file_to_override ${ARGN})
            get_filename_component(file_name ${original_file} NAME)
            if(${file_name} STREQUAL ${file_to_override})
                set(overridden_file
                    ${CMAKE_CURRENT_SOURCE_DIR}/${file_to_override})
            endif()
        endforeach()
        list(APPEND new_sources_list ${overridden_file})
    endforeach()
    set(${sources_list} ${new_sources_list} PARENT_SCOPE)
endfunction(override)


