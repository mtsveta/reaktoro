# Find Boost library
find_package(Boost REQUIRED)

# Find cvodes library
find_library(CVODES sundials_cvodes)
message("-- CVODES path: ${CVODES}")

# Find pybind11 library (if needed) https://github.com/pybind/pybind11
if(REAKTORO_BUILD_PYTHON)

    if(NOT PYBIND11_PYTHON_VERSION AND NOT PYTHON_EXECUTABLE)
        find_package(PythonInterp REQUIRED) # PYTHON_EXECUTABLE gets set here
    endif()

    find_package(pybind11 REQUIRED HINTS ${PROJECT_SOURCE_DIR}/thirdparty/pybind11)

    if(NOT pybind11_FOUND)
        message(WARNING "Could not find pybind11 - the Python module `reaktoro` will not be built.")
        set(REAKTORO_BUILD_PYTHON OFF)
    else()
        message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
    endif()

endif()

