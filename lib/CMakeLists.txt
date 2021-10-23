cmake_minimum_required(VERSION 3.4.3)
project(libiim CXX)

include_directories(${DIIM_SOURCE_DIR}/include)
include_directories(${Stdutils_INCLUDE_DIRS})
include_directories(${Numlib_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${DIIM_BINARY_DIR})
link_directories(${Numlib_LIBRARY_DIRS})
if(MKL_FOUND)
    link_directories(${MKL_LIBRARY_DIRS})
endif()

set(
    SRC_FILES
    diim.cpp
    perturbation.cpp
    utils.cpp
)
add_library(
    iim 
    ${SRC_FILES}
)
add_dependencies(iim Stdutils Numlib)
target_link_libraries(
    iim 
    ${BLAS_LIBRARIES}
    ${Numlib_LIBRARIES} 
) 

install(
    TARGETS iim
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
