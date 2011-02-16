set(SRC
    iso_c_utilities.f90
    cdirac.f90

    pyelk.cpp
)

find_package(Python REQUIRED)
find_package(NumPy REQUIRED)
find_package(Cython REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${PROJECT_BINARY_DIR}/src)

cython_add_module_pyx(pyelk)
add_python_library(pyelk ${SRC})
target_link_libraries(pyelk elkdirac)
