cmake_minimum_required(VERSION 3.18)
project(GPS_L1CA_GPU_SDR CUDA CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find CUDA
find_package(CUDA REQUIRED)
enable_language(CUDA)

# Set CUDA architecture (adjust based on your GPU)
set(CMAKE_CUDA_ARCHITECTURES 75)  # For RTX 2080, adjust as needed

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Library/C++_matrix_library
    ${CMAKE_CURRENT_SOURCE_DIR}/Library/C++_matrix_library/Eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/Library/fftw-3.3.5-dll64
    ${CUDA_INCLUDE_DIRS}
)

# Source files
set(SOURCES
    SDR_Development_CUDA.cpp
    AcqTrk.cpp
    DSP.cpp
    Global.cpp
    Init.cpp
    UI.cpp
    stdafx.cpp
)

# CUDA source files
set(CUDA_SOURCES
    CUDA_op.cu
)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES} ${CUDA_SOURCES})

# Link directories
link_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/Library/fftw-3.3.5-dll64
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link libraries
target_link_libraries(${PROJECT_NAME}
    ${CUDA_LIBRARIES}
    ${CUDA_cufft_LIBRARY}
    ${CUDA_cublas_LIBRARY}
    libfftw3-3
    libfftw3f-3
    libfftw3l-3
)

# Set CUDA properties
set_target_properties(${PROJECT_NAME} PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

# Copy DLLs to output directory (Windows)
if(WIN32)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/libfftw3-3.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/libfftw3f-3.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>)
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/libfftw3l-3.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>)
endif()

# Compiler flags
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W3)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)
endif()

# CUDA flags
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    --use_fast_math
    -Xcompiler "/wd4819"  # Disable warning about code page
>)