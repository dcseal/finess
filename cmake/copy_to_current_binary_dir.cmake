function(copy_to_current_binary_dir)
    file(COPY ${ARGN} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
