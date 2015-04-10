macro(check_python)
    find_package(PythonInterp REQUIRED)
    if(NOT ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} VERSION_EQUAL 2.7)
        message(WARNING "The Python scripts have only been tested with Python 2.7.")
    endif()
endmacro()

