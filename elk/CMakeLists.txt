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

set(PYTHON_INSTALL_PATH lib/python/site-packages)
set(PY_PATH ${PYTHON_INSTALL_PATH}/elk)
install(TARGETS pyelk
            RUNTIME DESTINATION ${PY_PATH}
            ARCHIVE DESTINATION ${PY_PATH}
            LIBRARY DESTINATION ${PY_PATH}
        )
install(FILES __init__.py mesh.py
    DESTINATION ${PY_PATH}
    )
