macro(setup_silo)
    find_package(SILO REQUIRED)
    if(SILO_FOUND)
        include_directories(${SILO_INCLUDE_DIR})
        list(APPEND EXTERNAL_LIBS ${SILO_LIBRARIES})
    endif()
    
    find_package(HDF5 REQUIRED)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIRS})
        list(APPEND EXTERNAL_LIBS ${HDF5_LIBRARIES})
    endif()
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS})
endmacro()
