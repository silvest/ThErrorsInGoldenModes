cmake_minimum_required(VERSION 3.15)
project(GoldenModes VERSION 1.0 LANGUAGES CXX)


# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find Boost
find_package(Boost REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost found: ${Boost_VERSION}")
endif()

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree MathCore Matrix)
if(ROOT_FOUND)
    message(STATUS "ROOT found: ${ROOT_VERSION}")
    include(${ROOT_USE_FILE})
endif()

find_package(MPI REQUIRED COMPONENTS CXX)

# Optional local BAT CMake build, keeping headers and binary ABI consistent.
set(BAT_SOURCE_DIR "" CACHE PATH "BAT source directory (for a local build)")
set(BAT_LIBRARY "" CACHE FILEPATH "BAT library from that local build")
if((BAT_SOURCE_DIR AND NOT BAT_LIBRARY) OR (BAT_LIBRARY AND NOT BAT_SOURCE_DIR))
    message(FATAL_ERROR "Set both BAT_SOURCE_DIR and BAT_LIBRARY for a local BAT build")
endif()
if(BAT_SOURCE_DIR AND BAT_LIBRARY)
    set(BAT_VERSION "local-PT")
    set(BAT_INCLUDE_DIRS "${BAT_SOURCE_DIR}")
    set(BAT_CFLAGS "")
    set(BAT_LIBS "${BAT_LIBRARY}")
    get_filename_component(BAT_RUNTIME_LIBRARY_DIR "${BAT_LIBRARY}" DIRECTORY)
else()
    find_package(BAT CONFIG QUIET)
    if(BAT_FOUND)
        set(BAT_INCLUDE_DIRS "${BAT_INCLUDE_DIR}")
        set(BAT_CFLAGS "")
        set(BAT_LIBS "BAT::BAT")
        set(BAT_RUNTIME_LIBRARY_DIR "${BAT_LIBRARY_DIR}")
        if(NOT BAT_HAS_PARALLEL_TEMPERING)
            message(FATAL_ERROR "The BAT package lacks Parallel Tempering; select the upgraded bat_MPI installation")
        endif()
    else()
        find_program(BAT_CONFIG_EXECUTABLE bat-config)
        if(NOT BAT_CONFIG_EXECUTABLE)
            message(FATAL_ERROR "bat-config not found; set BAT_SOURCE_DIR and BAT_LIBRARY for a local build")
        endif()
        execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --cflags OUTPUT_VARIABLE BAT_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --libs OUTPUT_VARIABLE BAT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --version OUTPUT_VARIABLE BAT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT BAT_VERSION)
            message(FATAL_ERROR "bat-config failed; check the BAT installation")
        endif()
        string(REGEX MATCHALL "-I[^ ]+" BAT_INCLUDE_FLAGS "${BAT_CFLAGS}")
        set(BAT_INCLUDE_DIRS "")
        execute_process(COMMAND ${BAT_CONFIG_EXECUTABLE} --libdir OUTPUT_VARIABLE BAT_RUNTIME_LIBRARY_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
        foreach(flag ${BAT_INCLUDE_FLAGS})
            string(SUBSTRING ${flag} 2 -1 include_dir)
            list(APPEND BAT_INCLUDE_DIRS "${include_dir}")
        endforeach()
    endif()
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${BAT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Source files for the library
set(LIBRARY_SOURCES
    goldenmodesB.cpp
    goldenmodesB_indSU3.cpp
    CKM.cpp
    ModelParameter.cpp
    PDGAverage.cpp
    CorrelatedGaussianObservables.cpp
    CorrelatedGaussianParameters.cpp
)

# Header files
set(LIBRARY_HEADERS
    goldenmodesB.h
    goldenmodesB_indSU3.h
    CKM.h
    ModelParameter.h
    PDGAverage.h
    CorrelatedGaussianObservables.h
    CorrelatedGaussianParameters.h
    dato.h
)

# Create main executable directly with all source files
add_executable(mainGM 
    mainGM.cpp
    goldenmodesB.cpp
    goldenmodesB_indSU3.cpp
    CKM.cpp
    ModelParameter.cpp
    CorrelatedGaussianObservables.cpp
    CorrelatedGaussianParameters.cpp
    PDGAverage.cpp
)

# Add compile options from BAT
separate_arguments(BAT_CFLAGS_LIST UNIX_COMMAND "${BAT_CFLAGS}")
target_compile_options(mainGM PRIVATE ${BAT_CFLAGS_LIST})

# Add BAT libraries
separate_arguments(BAT_LIBS_LIST UNIX_COMMAND "${BAT_LIBS}")

# Add link options to use libc++ (Clang/AppleClang only)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_options(mainGM PRIVATE -stdlib=libc++)
    if(APPLE)
        target_link_options(mainGM PRIVATE -L/opt/homebrew/opt/llvm/lib/c++)
    endif()
    target_link_libraries(mainGM 
        PRIVATE 
        ${BAT_LIBS_LIST}
        ${ROOT_LIBRARIES}
        c++
        c++abi
    )
else()
    target_link_libraries(mainGM 
        PRIVATE 
        ${BAT_LIBS_LIST}
        ${ROOT_LIBRARIES}
    )
endif()

if(TARGET MPI::MPI_CXX)
    target_link_libraries(mainGM PRIVATE MPI::MPI_CXX)
endif()

# Keep the selected BAT installation visible to the installed executable.
# Linux LD_LIBRARY_PATH can still override RUNPATH; avoid old BAT entries in jobs.
if(BAT_RUNTIME_LIBRARY_DIR)
    set_property(TARGET mainGM APPEND PROPERTY BUILD_RPATH "${BAT_RUNTIME_LIBRARY_DIR}")
    set_property(TARGET mainGM APPEND PROPERTY INSTALL_RPATH "${BAT_RUNTIME_LIBRARY_DIR}")
endif()

# Installation rules
install(TARGETS mainGM
    RUNTIME DESTINATION bin
)

install(FILES ${LIBRARY_HEADERS}
    DESTINATION include/GoldenModes
)

# Print configuration summary
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "  ROOT_VERSION: ${ROOT_VERSION}")
message(STATUS "  BAT_VERSION: ${BAT_VERSION}")
message(STATUS "  BAT_CFLAGS: ${BAT_CFLAGS}")
message(STATUS "  BAT_LIBS: ${BAT_LIBS}")
message(STATUS "")

option(GOLDENMODES_BUILD_TESTS "Build model prior-range regression test" OFF)
if(GOLDENMODES_BUILD_TESTS)
    enable_testing()
    add_executable(test_prior_ranges tests/prior_ranges.cpp ${LIBRARY_SOURCES})
    target_compile_options(test_prior_ranges PRIVATE ${BAT_CFLAGS_LIST})
    target_link_libraries(test_prior_ranges PRIVATE ${BAT_LIBS_LIST} ${ROOT_LIBRARIES})
    if(TARGET MPI::MPI_CXX)
        target_link_libraries(test_prior_ranges PRIVATE MPI::MPI_CXX)
    endif()
    add_test(NAME prior_ranges COMMAND test_prior_ranges)
endif()
