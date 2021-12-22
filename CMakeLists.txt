cmake_minimum_required(VERSION 3.13)
project(DIIM 
    VERSION 0.0.1
    LANGUAGES CXX
)

################################################################################

include(GNUInstallDirs)

################################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

################################################################################

set(Diim_STANDALONE_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(Diim_STANDALONE_PROJECT ON)
endif()

option(Diim_BUILD_TESTS "Build tests." ${Diim_STANDALONE_PROJECT})

################################################################################

# Set installation directory.
if (WIN32)
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX $ENV{USERPROFILE} CACHE PATH "USERPROFILE" FORCE)
    endif()
else()
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE PATH "HOME" FORCE)
    endif()
endif()

################################################################################

if(WIN32)
    option(BUILD_SHARED_LIBS "Compile shared libraries." OFF)
else()
    option(BUILD_SHARED_LIBS "Compile shared libraries." ON)
endif()

option(Diim_ENABLE_CONCEPTS "Try to enable concepts support by giving extra flags." ON)

# Check if C++20 is available.
if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(STATUS "Using C++20 standard")
    set(CMAKE_CXX_STANDARD 20)
else()
    message(FATAL_ERROR "Requested DIIM_CXX_STANDARD \"20\" not supported by provided C++ compiler")
endif()

if(Diim_ENABLE_CONCEPTS)
    if(CMAKE_CXX_STANDARD STREQUAL "20")
        include(CheckCXXCompilerFlag)
        CHECK_CXX_COMPILER_FLAG("-fconcepts" COMPILER_SUPPORTS_FCONCEPTS)
        if(COMPILER_SUPPORTS_FCONCEPTS)
            message(STATUS "-- Using \"-fconcepts\" to enable concepts support")
            set(Diim_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fconcepts)
        else()
            CHECK_CXX_COMPILER_FLAG("-fconcepts-ts" COMPILER_SUPPORTS_FCONCEPTS_TS)
            if(COMPILER_SUPPORTS_FCONCEPTS)
                message(STATUS "-- Using \"-fconcepts-ts\" to enable concepts support")
                set(Diim_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fconcepts-ts)
            endif()
        endif()
    endif()
endif()

# Set compiler options.
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(Diim_CXX_FLAGS ${Diim_CXX_FLAGS} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow -Wno-error=unused-parameter -Wno-maybe-uninitialized)
    set(Diim_CXX_FLAGS_DEBUG ${Diim_CXX_FLAGS} -g)
    set(Diim_CXX_FLAGS_RELEASE ${Diim_CXX_FLAGS} -O3)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(Diim_CMAKE_CXX_FLAGS ${Diim_CXX_FLAGS} -Wpedantic -Wall -Wshadow -Wextra -Wno-missing-braces -Werror -Wno-error=strict-overflow)
    set(Diim_CXX_FLAGS_DEBUG ${Diim_CXX_FLAGS} -g)
    set(Diim_CXX_FLAGS_RELEASE ${Diim_CXX_FLAGS} -O3)
endif()
if(MSVC)
    if (CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        set(Diim_CXX_FLAGS ${Diim_CXX_FLAGS} /W4 /Qstd:c++20)
        set(Diim_CXX_FLAGS_DEBUG ${Diim_CXX_FLAGS} /Zi /Od)
        set(Diim_CXX_FLAGS_RELEASE ${Diim_CXX_FLAGS} /O3)
    else()
        set(Diim_CXX_FLAGS ${Diim_CXX_FLAGS} /W4 /WX /wd4100)
        set(Diim_CXX_FLAGS_DEBUG ${Diim_CXX_FLAGS} /Zi /Od)
        set(Diim_CXX_FLAGS_RELEASE ${Diim_CXX_FLAGS} /O2)
    endif()
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

################################################################################

# mdspan is required.
find_package(mdspan QUIET)
if (NOT mdspan_FOUND)
    message(STATUS "No installed mdspan found, fetching from Github")
    include(FetchContent)
    FetchContent_Declare(
        mdspan
        GIT_REPOSITORY  https://github.com/kokkos/mdspan.git
        GIT_TAG         stable
    )
    FetchContent_MakeAvailable(mdspan)
endif()

if(NOT TARGET std::mdspan)
    find_package(mdspan REQUIRED)
endif()

# Scilib is required.
find_package(scilib QUIET)
if(scilib_FOUND)
    message(STATUS "Scilib found")
endif()
if (NOT scilib_FOUND)
    message(STATUS "No installed Scilib found, fetching from Github")
    include(FetchContent)
    FetchContent_Declare(
        scilib
        GIT_REPOSITORY https://github.com/stigrs/scilib.git
        GIT_TAG        main
    )
    FetchContent_MakeAvailable(scilib)
endif()

if(NOT TARGET scilib::scilib)
  find_package(scilib REQUIRED)
endif()

# Stdutils is required.
find_package(stdutils QUIET)
if(Diim_FOUND)
    message(STATUS "Stdutils found")
endif()
if (NOT Diim_FOUND)
    message(STATUS "No installed Stdutils found, fetching from Github")
    include(FetchContent)
    FetchContent_Declare(
        stdutils
        GIT_REPOSITORY https://github.com/stigrs/stdutils.git
        GIT_TAG        master
    )
    FetchContent_MakeAvailable(stdutils)
endif()

if(NOT TARGET stdutils::stdutils)
  find_package(stdutils REQUIRED)
endif()

# BLAS and LAPACKE are required.
if(DEFINED ENV{MKLROOT})
	find_package(MKL)
	if(MKL_FOUND)
	    message(STATUS "Intel MKL found")
		add_definitions(${MKL_DEFINITIONS})
	endif()
else()
    find_package(OpenBLAS)
	if(BLAS_FOUND)
        MESSAGE(STATUS "OpenBLAS found")
		add_definitions(${BLAS_DEFINITIONS})
	endif()
endif()

################################################################################

set(SRC_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/diim.cpp    
    ${CMAKE_CURRENT_SOURCE_DIR}/src/perturbation.cpp    
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp    
)

add_library(diim ${SRC_FILES})
add_library(diim::diim ALIAS diim)

target_link_libraries(
    diim 
    PRIVATE
    scilib::scilib
    stdutils::stdutils
    std::mdspan
    ${BLAS_LIBRARIES}
)

target_compile_options(
    diim
    PRIVATE
    "$<$<CONFIG:Debug>:${Diim_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${Diim_CXX_FLAGS_RELEASE}>"
)

target_include_directories(diim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${BLAS_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${LAPACKE_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS diim
    EXPORT diimTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}    
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}    
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}    
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}    
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/diim DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT diimTargets
    FILE diimTargets.cmake
    NAMESPACE diim::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/diim
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/diimConfigVersion.cmake"
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/diimConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/diimConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/diim
)
install(FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/diimConfig.cmake" 
    "${CMAKE_CURRENT_BINARY_DIR}/diimConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/diim
)

################################################################################

add_subdirectory(tools)

################################################################################

if(Diim_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
