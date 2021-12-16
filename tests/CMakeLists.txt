cmake_minimum_required(VERSION 3.17)
project(DIIM_tests CXX)

include_directories(${BLAS_INCLUDE_DIRS})
include_directories(${LAPACKE_INCLUDE_DIRS})
link_directories(${BLAS_LIBRARY_DIRS})

# Copy *.inp files to build directory.
file(GLOB INP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.inp")
foreach(inp_file_name ${INP_FILES})
    file(COPY "${inp_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

# Copy *.csv files to build directory.
file(GLOB CSV_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.csv")
foreach(csv_file_name ${CSV_FILES})
    file(COPY "${csv_file_name}" DESTINATION "${CMAKE_BINARY_DIR}/tests")
endforeach()

include(FetchContent)
FetchContent_Declare(
    googletest
    # Specify the commit you depend on and update it regularly.
    URL https://github.com/google/googletest/archive/refs/heads/main.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

if(NOT TARGET std::mdspan)
    find_package(mdspan REQUIRED)
endif()

if(NOT TARGET Scilib::Scilib)
    find_package(scilib REQUIRED)
endif()

if(NOT TARGET Stdutils::Stdutils)
    find_package(stdutils REQUIRED)
endif()

if(NOT TARGET Diim::Diim)
    find_package(diim REQUIRED)
endif()

function(add_diim_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries (
        ${name} 
        PRIVATE
        Diim::Diim
        Scilib::Scilib
        Stdutils::Stdutils
        gtest_main
        std::mdspan
	      ${BLAS_LIBRARIES}
    ) 
    add_test(NAME ${name} COMMAND ${name})
endfunction()

set(PROGRAMS 
    test_diim
)

foreach(program ${PROGRAMS})
    add_diim_test(${program})
	  if(WIN32)
	      if(NOT MKL_FOUND)
            set_tests_properties(${program} PROPERTIES ENVIRONMENT "PATH=${BLAS_BINARY_DIRS};$ENV{PATH}")
	      endif()
    else()
	      if(NOT MKL_FOUND)
            set_tests_properties(${program} PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:${BLAS_BINARY_DIRS}")
	      endif()
	  endif()
endforeach()
