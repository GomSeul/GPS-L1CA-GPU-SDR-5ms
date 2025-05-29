# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a GPS L1 C/A Software Defined Radio (SDR) implementation with GPU acceleration using CUDA. The project processes GPS satellite signals to compute position, velocity, and time (PVT) solutions with carrier phase measurements.

## Build Instructions

**Note:** This appears to be a Visual Studio project for Windows, but the solution files are not present in the repository. The project requires:

1. **NVIDIA CUDA Toolkit** - For GPU acceleration
2. **Visual Studio** with C++ support
3. **Dependencies** (already included):
   - FFTW 3.3.5 (DLLs included)
   - Eigen (source included in Library/C++_matrix_library/)

To build:
- Create a Visual Studio project including all .cpp, .cu, .h, and .cuh files
- Link against CUDA libraries (cudart, cufft, cublas)
- Link against FFTW libraries (libfftw3-3.dll included)
- Set up include paths for Eigen library

## Architecture

### Core Components

1. **SDR_Development_CUDA.cpp** - Main entry point, orchestrates the SDR processing pipeline
2. **CUDA_op.cu/cuh** - GPU-accelerated operations for acquisition and tracking
3. **AcqTrk.cpp/h** - Acquisition and tracking algorithms
4. **DSP.cpp/h** - Digital signal processing functions
5. **Init.cpp/h** - System initialization
6. **Global.cpp/h** - Global variables and shared state
7. **UI.cpp/h** - User interface and display functions

### Key Processing Parameters (from Defines.h)

- Sampling frequency: 25 MHz
- IF frequency: 5.42 MHz  
- Thread configuration: 512 threads per block
- Acquisition threshold: 5.5
- Tracking loops:
  - FLL bandwidth: 10 Hz
  - PLL bandwidth: 5 Hz
  - DLL bandwidth: 0.25 Hz

### Data Flow

1. Read IF data samples
2. Parallel acquisition on GPU to find satellites
3. Track satellites using FLL-assisted PLL and DLL
4. Decode navigation data
5. Compute PVT solution with corrections
6. Output results to DataLog/ directory

## Testing and Validation

MATLAB scripts in `DataLog/MATLAB_Plot/` are used for:
- Plotting acquisition results
- Visualizing tracking performance
- Comparing navigation solutions
- Analyzing trajectories

Test data files:
- `1_Motion_scenario_100_sec_stop.bin`
- `2_Motion_scenario_80_stop_and_turn.bin`
- `3_Motion_scenario_complex.bin`

## Output Files

Navigation results are saved to:
- `DataLog/DataLog_Navi_Out.txt` - ECEF coordinates
- `DataLog/DataLog_Navi_LLH_Out.txt` - Lat/Lon/Height coordinates
- GPU versions with `_GPU` suffix when GPU acceleration is used

## Reference Documentation

The project includes the following reference material in the `doc/` folder:

**"Understanding GPS/GNSS: Principles and Applications (3rd Edition)"** by Elliott D. Kaplan and Christopher J. Hegarty
- This comprehensive textbook should be consulted for:
  - GPS signal structure and modulation details
  - Acquisition and tracking loop theory
  - Navigation solution algorithms
  - Error sources and correction methods
  - Carrier phase measurements and ambiguity resolution

When implementing or modifying GPS/GNSS algorithms in this codebase, refer to this textbook for theoretical background and implementation guidance.