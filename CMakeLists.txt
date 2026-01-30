cmake_minimum_required(VERSION 3.15)
project(GoldenModes VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_COMPILER "mpicxx")

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree MathCore Matrix)
if(ROOT_FOUND)
    message(STATUS "ROOT found: ${ROOT_VERSION}")
    include(${ROOT_USE_FILE})
endif()

# Find BAT (Bayesian Analysis Toolkit) using bat-config
execute_process(
    COMMAND bat-config --cflags
    OUTPUT_VARIABLE BAT_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND bat-config --libs
    OUTPUT_VARIABLE BAT_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND bat-config --version
    OUTPUT_VARIABLE BAT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(BAT_VERSION)
    message(STATUS "BAT found: ${BAT_VERSION}")
    # Extract include directories from CFLAGS
    string(REGEX MATCHALL "-I[^ ]+" BAT_INCLUDE_FLAGS ${BAT_CFLAGS})
    set(BAT_INCLUDE_DIRS "")
    foreach(flag ${BAT_INCLUDE_FLAGS})
        string(SUBSTRING ${flag} 2 -1 include_dir)
        list(APPEND BAT_INCLUDE_DIRS ${include_dir})
    endforeach()
else()
    message(FATAL_ERROR "BAT not found. Please install BAT or ensure bat-config is in PATH.")
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${BAT_INCLUDE_DIRS}
)

# Source files for the library
set(LIBRARY_SOURCES
    goldenmodesB.cpp
    CKM.cpp
    ModelParameter.cpp
    PDGAverage.cpp
    CorrelatedGaussianObservables.cpp
    CorrelatedGaussianParameters.cpp
)

# Header files
set(LIBRARY_HEADERS
    goldenmodesB.h
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

# Add link options to use libc++
target_link_options(mainGM PRIVATE -stdlib=libc++ -L/opt/homebrew/opt/llvm/lib/c++)

target_link_libraries(mainGM 
    PRIVATE 
    ${BAT_LIBS_LIST}
    ${ROOT_LIBRARIES}
    c++
    c++abi
)

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
