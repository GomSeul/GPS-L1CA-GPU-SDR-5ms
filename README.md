# GPS L1 C/A GPU-Accelerated Software Defined Radio with Carrier Phase

A GPU-accelerated GPS L1 C/A signal processing software defined radio (SDR) implementation with carrier phase measurements for precise positioning.

## Features

- GPU-accelerated signal acquisition and tracking using NVIDIA CUDA
- GPS L1 C/A signal processing with carrier phase measurements
- Real-time signal processing at 25 MHz sampling rate
- 5.42 MHz intermediate frequency (IF)
- Parallel acquisition across multiple satellites
- FLL-assisted PLL tracking (FLL: 10 Hz, PLL: 5 Hz, DLL: 0.25 Hz)
- Support for both 1ms and 5ms coherent integration
- Navigation solution with corrections:
  - Earth rotation
  - Satellite clock
  - Ionospheric delay
  - Tropospheric delay

## Requirements

- Windows OS (Visual Studio project)
- NVIDIA GPU with CUDA support
- CUDA Toolkit
- Visual Studio with C++ support
- Dependencies (included):
  - FFTW 3.3.5 (Fast Fourier Transform library)
  - Eigen (C++ matrix library)

## Project Structure

```
├── CUDA_op.cu/cuh      # GPU kernels for acquisition and tracking
├── AcqTrk.cpp/h        # Acquisition and tracking algorithms
├── DSP.cpp/h           # Digital signal processing functions
├── Init.cpp/h          # System initialization
├── Global.cpp/h        # Global variables and shared state
├── UI.cpp/h            # User interface
├── Defines.h           # Configuration parameters
├── Structs.h           # Data structures
├── CLAUDE.md           # Development guidelines
├── DataLog/            # Output navigation data
│   └── MATLAB_Plot/    # MATLAB visualization scripts
├── Library/            # External libraries
│   └── C++_matrix_library/  # Eigen library
└── doc/                # Reference documentation
```

## Key Parameters

- Sampling frequency: 25 MHz
- IF frequency: 5.42 MHz
- Thread configuration: 512 threads per block
- Acquisition threshold: 5.5
- Tracking threshold: 1500

## Recent Updates

- Enhanced acquisition to support 5ms coherent integration
- Added CUDA kernels for 5ms processing:
  - `f_Replica_Gen_5ms()` - Generate 5ms replicas
  - `f_Freq_Conjugate_5ms()` - Frequency domain conjugate multiplication
  - `f_Kernel_Abs_5ms()` - Calculate absolute values
  - `f_Acquisition_5ms()` - Main 5ms acquisition function
  - `f_Peak_Detector_5ms()` - Peak detection for 5ms data

## Build Instructions

1. Open project in Visual Studio
2. Ensure CUDA Toolkit is installed and configured
3. Add all .cpp, .cu, .h, and .cuh files to the project
4. Configure project properties:
   - Set CUDA include directories
   - Link against CUDA libraries (cudart, cufft, cublas)
   - Link against FFTW libraries (included DLLs)
   - Add Eigen include path
5. Build the solution

## Output Files

The SDR outputs navigation results to:
- `DataLog/DataLog_Navi_Out.txt` - ECEF coordinates
- `DataLog/DataLog_Navi_LLH_Out.txt` - Latitude/Longitude/Height
- GPU-accelerated versions saved with `_GPU` suffix

## Testing

Use the included MATLAB scripts in `DataLog/MATLAB_Plot/` to visualize:
- Acquisition results
- Tracking performance
- Navigation trajectories
- Velocity profiles

Test data files included:
- `1_Motion_scenario_100_sec_stop.bin`
- `2_Motion_scenario_80_stop_and_turn.bin`
- `3_Motion_scenario_complex.bin`

## Reference

For theoretical background and implementation details, refer to:
- "Understanding GPS/GNSS: Principles and Applications (3rd Edition)" by Elliott D. Kaplan and Christopher J. Hegarty (included in doc/ folder)

## License

This project is for educational and research purposes.

## Contributors

- Original implementation with CUDA GPU acceleration
- Enhanced with 5ms coherent integration support