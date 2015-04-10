macro(check_compiler)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.8.2)
            message(WARNING "This program has only been tested with gcc 4.8.2.")
        endif()
        set(USE_GCC YES)
    else()
        message(WARNING "This program has only been tested with gcc 4.8.2, and may contain code incompatible with non-gcc compilers.")
    endif()
endmacro()
