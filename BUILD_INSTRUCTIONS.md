# Build Instructions for GPS L1CA GPU SDR

## Prerequisites

1. **Visual Studio 2019 or later** with C++ support
2. **NVIDIA CUDA Toolkit** (version 11.0 or later recommended)
3. **NVIDIA GPU** with compute capability 3.5 or higher

## Build Steps in Visual Studio

### 1. Create New Project
1. Open Visual Studio
2. Create New Project → CUDA Runtime Project
3. Name: GPS_L1CA_GPU_SDR
4. Location: Choose your preferred directory

### 2. Add Existing Files
1. Right-click on project → Add → Existing Item
2. Add all source files:
   - `SDR_Development_CUDA.cpp`
   - `AcqTrk.cpp`
   - `DSP.cpp`
   - `Global.cpp`
   - `Init.cpp`
   - `UI.cpp`
   - `stdafx.cpp`
   - `CUDA_op.cu`

3. Add all header files:
   - `AcqTrk.h`
   - `CUDA_op.cuh`
   - `DSP.h`
   - `Defines.h`
   - `Global.h`
   - `Init.h`
   - `Structs.h`
   - `UI.h`
   - `stdafx.h`
   - `targetver.h`

### 3. Configure Project Properties

Right-click on project → Properties:

#### General
- Configuration Type: Application (.exe)
- Platform Toolset: Latest available

#### C/C++ → General
- Additional Include Directories:
  ```
  $(ProjectDir)Library\C++_matrix_library
  $(ProjectDir)Library\C++_matrix_library\Eigen
  $(ProjectDir)Library\fftw-3.3.5-dll64
  $(CudaToolkitIncludeDir)
  ```

#### C/C++ → Preprocessor
- Preprocessor Definitions:
  ```
  WIN32
  _DEBUG (for Debug) or NDEBUG (for Release)
  _CONSOLE
  _UNICODE
  UNICODE
  ```

#### CUDA C/C++ → Common
- Additional Include Directories: Same as C/C++ includes
- Generate Relocatable Device Code: Yes (-rdc=true)

#### CUDA C/C++ → Device
- Code Generation: compute_75,sm_75 (adjust based on your GPU)

#### Linker → General
- Additional Library Directories:
  ```
  $(ProjectDir)Library\fftw-3.3.5-dll64
  $(CudaToolkitLibDir)
  ```

#### Linker → Input
- Additional Dependencies:
  ```
  cudart.lib
  cufft.lib
  cublas.lib
  libfftw3-3.lib
  libfftw3f-3.lib
  libfftw3l-3.lib
  ```

### 4. Copy DLL Files
Copy the following DLLs to your output directory (Debug/Release):
- `libfftw3-3.dll`
- `libfftw3f-3.dll`
- `libfftw3l-3.dll`

### 5. Build
1. Select configuration (Debug/Release) and platform (x64)
2. Build → Build Solution (F7)

## Troubleshooting

### Common Issues

1. **CUDA not found**
   - Ensure CUDA Toolkit is installed
   - Check environment variable: `CUDA_PATH`

2. **Cannot open include file**
   - Verify all include directories are correctly set
   - Check that Eigen library path is correct

3. **Unresolved external symbols**
   - Ensure all CUDA libraries are linked
   - Check that FFTW libraries are in the correct path

4. **DLL not found**
   - Copy FFTW DLLs to the same directory as the executable

### GPU Architecture
Adjust the code generation settings based on your GPU:
- GTX 1060/1070/1080: compute_61,sm_61
- RTX 2060/2070/2080: compute_75,sm_75
- RTX 3060/3070/3080/3090: compute_86,sm_86
- RTX 4060/4070/4080/4090: compute_89,sm_89

## Running the Application

1. Ensure you have GPS IF data file
2. Run the executable from command prompt or Visual Studio
3. Monitor the output in DataLog directory