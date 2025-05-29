#pragma once
#include <cuda_runtime.h>
#include "cuda.h"
#include <cufft.h>
//#include <cufftw.h>
#include "cublas_v2.h"
#include "Defines.h"
#include "device_launch_parameters.h"



void f_Acquisition(
	short* PRN,						// ���� PRN
	short* ifdata,					// if ������
	int* CA_code_GPS_1023,
	short* dopp_freq,
	cufftDoubleComplex* replica_carr_vec,
	cufftDoubleComplex* replica_code_vec,
	cufftDoubleComplex* fft_out_carr_vec,
	cufftDoubleComplex* fft_out_code_vec,
	cufftDoubleComplex* freq_conj_carr_code,
	cufftDoubleComplex* dump_IFFT,
	double* Z_IFFT,
	double* Cpu_Z_IFFT,
	int* ACQ_IF_DUMP,
	int* ACQ_MODE_V2,
	int* ACTIVECHANNEL,
	int* Cpu_code_phase
);

void f_Acquisition_5ms(
	short* PRN,						// ���� PRN
	short* ifdata,					// if ������ (5ms)
	int* CA_code_GPS_1023,
	short* dopp_freq,
	cufftDoubleComplex* replica_carr_vec,
	cufftDoubleComplex* replica_code_vec,
	cufftDoubleComplex* fft_out_carr_vec,
	cufftDoubleComplex* fft_out_code_vec,
	cufftDoubleComplex* freq_conj_carr_code,
	cufftDoubleComplex* dump_IFFT,
	double* Z_IFFT,
	double* Cpu_Z_IFFT,
	int* ACQ_IF_DUMP,
	int* ACQ_MODE_V2,
	int* ACTIVECHANNEL,
	int* Cpu_code_phase
);

void f_constant_setup(int setup_f_s, int setup_f_if_GPS, int setup_f_ca_GPS, unsigned int sample_size_1ms);

void f_Correlation
(
	short* PRN,						// ���� PRN
	short* code_index,
	short* carr_temp,
	int* s_c,
	short* ifdata,					// if ������
	float* dump,
	int* cuda_CA_code_GPS_2046,
	int ACTIVECHANNEL,
	int* trk_CH_2046,
	int TRK_first,
	int* trk_count
);

