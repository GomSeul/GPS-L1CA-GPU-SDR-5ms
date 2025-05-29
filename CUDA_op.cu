// CUDA compilation fix
#ifdef __CUDACC__
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include "CUDA_op.cuh"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include "Global.h"

// Forward declarations to avoid including AcqTrk.h
void f_Peak_Detector(Channel_struct& f_CH, double* Z_IFFT);
void f_Peak_Detector_5ms(Channel_struct& f_CH, double* Z_IFFT);

// #include "AcqTrk.h"  // Commented out to avoid stdafx.h dependency


#define CorrSizeLocal SAMPLING_FREQ / 1000.0		//�⺻�� : 50000
#define Trk_CorrSizeLocal (SAMPLING_FREQ / 1000.0) + 10		//�⺻�� : 50000
#define BATCH 1
#define SAMPLING_FREQ	Sampling_Frequency
#define DIM_BLOCK_COL	512
#define F_CA_GPS		1023000
#define IF_FREQ			IF_Frequency

__constant__ int Gpu_const_f_s[1] = { 0, };
__constant__ int Gpu_const_f_if_GPS[1] = { 0, };
__constant__ int Gpu_const_f_ca_GPS[1] = { 0, };
__constant__ int Gpu_const_sample_size_1ms[1] = { 0, };


__global__ void f_Replica_Gen(cufftDoubleComplex* replica_carr_vec, cufftDoubleComplex* replica_code_vec, short* IF_Data, short* Range_dopp, int* CA_code_GPS_1023, short* PRN)
{
	int thread_num = blockDim.x * blockIdx.x + threadIdx.x;
	int CH_PRN = blockIdx.y / 45;
	short ifdata_shift = IF_Data[thread_num] >> 8;

	if (thread_num > (int)(SAMPLING_FREQ / 1000.0))	return;

	int dopp_f = Range_dopp[blockIdx.y % 45];

	int temp = ((int)((double)thread_num / (double)SAMPLING_FREQ * (double)F_CA_GPS)) % 1023;
	int element_num = thread_num + ((int)(SAMPLING_FREQ / 1000.0) * blockIdx.y);

	replica_carr_vec[element_num].x = ifdata_shift * cos(2.0 * M_PI * (IF_FREQ + dopp_f) * thread_num / SAMPLING_FREQ);
	replica_carr_vec[element_num].y = ifdata_shift * sin(2.0 * M_PI * (IF_FREQ + dopp_f) * thread_num / SAMPLING_FREQ);
	replica_code_vec[element_num].x = CA_code_GPS_1023[(PRN[CH_PRN] - 1) * 1023 + temp];
	replica_code_vec[element_num].y = 0.0;
}

__global__ void f_Replica_Gen_5ms(cufftDoubleComplex* replica_carr_vec, cufftDoubleComplex* replica_code_vec, short* IF_Data, short* Range_dopp, int* CA_code_GPS_1023, short* PRN)
{
	int thread_num = blockDim.x * blockIdx.x + threadIdx.x;
	int CH_PRN = blockIdx.y / 45;
	
	if (thread_num >= (int)(SAMPLING_FREQ * 5 / 1000.0))	return;  // 5ms boundary check
	
	short ifdata_shift = IF_Data[thread_num] >> 8;
	int dopp_f = Range_dopp[blockIdx.y % 45];

	// For 5ms, we need to repeat the CA code 5 times (5 x 1023 chips)
	int temp = ((int)((double)thread_num / (double)SAMPLING_FREQ * (double)F_CA_GPS)) % 1023;
	int element_num = thread_num + ((int)(SAMPLING_FREQ * 5 / 1000.0) * blockIdx.y);

	replica_carr_vec[element_num].x = ifdata_shift * cos(2.0 * M_PI * (IF_FREQ + dopp_f) * thread_num / SAMPLING_FREQ);
	replica_carr_vec[element_num].y = ifdata_shift * sin(2.0 * M_PI * (IF_FREQ + dopp_f) * thread_num / SAMPLING_FREQ);
	replica_code_vec[element_num].x = CA_code_GPS_1023[(PRN[CH_PRN] - 1) * 1023 + temp];
	replica_code_vec[element_num].y = 0.0;
}

__global__ void f_Freq_Conjugate(cufftDoubleComplex* freq_Mul_conj_carr_code, cufftDoubleComplex* freq_replica_carr_vec, cufftDoubleComplex* freq_replica_code_vec)
{
	int thread_num = blockDim.x * blockIdx.x + threadIdx.x;
	if (thread_num > (int)(SAMPLING_FREQ / 1000.0))	return;

	int element_num = thread_num + ((int)(SAMPLING_FREQ / 1000.0) * blockIdx.y);

	freq_Mul_conj_carr_code[element_num].x = (freq_replica_code_vec[element_num].x * freq_replica_carr_vec[element_num].x) + (freq_replica_code_vec[element_num].y * freq_replica_carr_vec[element_num].y);
	freq_Mul_conj_carr_code[element_num].y = (freq_replica_code_vec[element_num].x * freq_replica_carr_vec[element_num].y) - (freq_replica_code_vec[element_num].y * freq_replica_carr_vec[element_num].x);


}

__global__ void f_Freq_Conjugate_5ms(cufftDoubleComplex* freq_Mul_conj_carr_code, cufftDoubleComplex* freq_replica_carr_vec, cufftDoubleComplex* freq_replica_code_vec)
{
	int thread_num = blockDim.x * blockIdx.x + threadIdx.x;
	if (thread_num >= (int)(SAMPLING_FREQ * 5 / 1000.0))	return;  // 5ms boundary check

	int element_num = thread_num + ((int)(SAMPLING_FREQ * 5 / 1000.0) * blockIdx.y);

	freq_Mul_conj_carr_code[element_num].x = (freq_replica_code_vec[element_num].x * freq_replica_carr_vec[element_num].x) + (freq_replica_code_vec[element_num].y * freq_replica_carr_vec[element_num].y);
	freq_Mul_conj_carr_code[element_num].y = (freq_replica_code_vec[element_num].x * freq_replica_carr_vec[element_num].y) - (freq_replica_code_vec[element_num].y * freq_replica_carr_vec[element_num].x);
}


__global__ void f_Kernel_Abs(cufftDoubleComplex* dump_IFFT, double* Z_IFFT)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int bch = blockIdx.y;	//doppler frequency
	int index = (bch * CorrSizeLocal) + tid;

	while (tid < (int)(*Gpu_const_sample_size_1ms))
	{
		Z_IFFT[index] = (dump_IFFT[index].x * dump_IFFT[index].x) + (dump_IFFT[index].y * dump_IFFT[index].y);

		tid += blockDim.x * gridDim.x;
	}
	__syncthreads();

}

__global__ void f_Kernel_Abs_5ms(cufftDoubleComplex* dump_IFFT, double* Z_IFFT)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int bch = blockIdx.y;	//doppler frequency
	int sample_size_5ms = (int)(SAMPLING_FREQ * 5 / 1000.0);
	int index = (bch * sample_size_5ms) + tid;

	while (tid < sample_size_5ms)
	{
		Z_IFFT[index] = (dump_IFFT[index].x * dump_IFFT[index].x) + (dump_IFFT[index].y * dump_IFFT[index].y);

		tid += blockDim.x * gridDim.x;
	}
	__syncthreads();
}

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
)
{

	int c = 0;

	dim3 grids(blocksPerGrid, doppler_bin_size * 32, 1);

	cufftHandle p_cufft1, p_cufft2, p_cuifft;

	f_Replica_Gen << <grids, threadsPerBlock >> > (replica_carr_vec, replica_code_vec, ifdata, dopp_freq, CA_code_GPS_1023, PRN);

	int sampling_freq = SAMPLING_FREQ;
	int n_doppler_bin = doppler_bin_size;
	int batch = doppler_bin_size * 32;


	//FFT Plan
	cufftPlan1d(&p_cufft1, (int)(SAMPLING_FREQ / 1000.0), CUFFT_Z2Z, batch);
	cufftPlan1d(&p_cufft2, (int)(SAMPLING_FREQ / 1000.0), CUFFT_Z2Z, batch);
	cufftPlan1d(&p_cuifft, (int)(SAMPLING_FREQ / 1000.0), CUFFT_Z2Z, batch);


	//FFT ����
	cufftExecZ2Z(p_cufft1, replica_carr_vec, fft_out_carr_vec, CUFFT_FORWARD);
	cufftExecZ2Z(p_cufft2, replica_code_vec, fft_out_code_vec, CUFFT_FORWARD);

	cudaDeviceSynchronize();

	f_Freq_Conjugate << <grids, threadsPerBlock >> > (freq_conj_carr_code, fft_out_carr_vec, fft_out_code_vec);

	//IFFT Plan
	cufftExecZ2Z(p_cuifft, freq_conj_carr_code, dump_IFFT, CUFFT_INVERSE);

	cudaDeviceSynchronize();


	f_Kernel_Abs << <grids, threadsPerBlock >> > (dump_IFFT, Z_IFFT);


	cudaMemcpy(Cpu_Z_IFFT, Z_IFFT, sizeof(double) * CorrSizeLocal * doppler_bin_size * 32, cudaMemcpyDeviceToHost);

	while ((*ACQ_IF_DUMP) == 1 && (*ACQ_MODE_V2) == 1)
	{

		f_Peak_Detector(CH[c], Cpu_Z_IFFT);

		if (CH[c].detect == true)
		{
			(*ACTIVECHANNEL)++;
			Cpu_code_phase[c] = CH[c].codephase;

			if (start_mode == false)
			{
				if (CH[c].PRN < 32)
				{
					CH[c + 1].PRN = CH[c].PRN + 1;
					c++;

				}
				else
				{
					(*ACQ_MODE_V2) = 0;
				}
			}
			else
			{
				c++;
			}
		}
		if (CH[c].ACQ_Check == false)
		{
			(*ACQ_MODE_V2) = 0;

		}
	}


	cufftDestroy(p_cufft1);
	cufftDestroy(p_cufft2);
	cufftDestroy(p_cuifft);


	//return true;
}

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
)
{
	int c = 0;
	int sample_size_5ms = (int)(SAMPLING_FREQ * 5 / 1000.0);

	// Use 5ms grid configuration
	dim3 grids(blocksPerGrid_5ms, doppler_bin_size * 32, 1);

	cufftHandle p_cufft1, p_cufft2, p_cuifft;

	// Generate 5ms replicas
	f_Replica_Gen_5ms<<<grids, threadsPerBlock>>>(replica_carr_vec, replica_code_vec, ifdata, dopp_freq, CA_code_GPS_1023, PRN);

	int batch = doppler_bin_size * 32;

	// FFT Plan for 5ms data
	cufftPlan1d(&p_cufft1, sample_size_5ms, CUFFT_Z2Z, batch);
	cufftPlan1d(&p_cufft2, sample_size_5ms, CUFFT_Z2Z, batch);
	cufftPlan1d(&p_cuifft, sample_size_5ms, CUFFT_Z2Z, batch);

	// FFT execution
	cufftExecZ2Z(p_cufft1, replica_carr_vec, fft_out_carr_vec, CUFFT_FORWARD);
	cufftExecZ2Z(p_cufft2, replica_code_vec, fft_out_code_vec, CUFFT_FORWARD);

	cudaDeviceSynchronize();

	// Frequency domain conjugate multiplication for 5ms
	f_Freq_Conjugate_5ms<<<grids, threadsPerBlock>>>(freq_conj_carr_code, fft_out_carr_vec, fft_out_code_vec);

	// IFFT
	cufftExecZ2Z(p_cuifft, freq_conj_carr_code, dump_IFFT, CUFFT_INVERSE);

	cudaDeviceSynchronize();

	// Calculate absolute values for 5ms
	f_Kernel_Abs_5ms<<<grids, threadsPerBlock>>>(dump_IFFT, Z_IFFT);

	// Copy results back to host (5ms size)
	cudaMemcpy(Cpu_Z_IFFT, Z_IFFT, sizeof(double) * sample_size_5ms * doppler_bin_size * 32, cudaMemcpyDeviceToHost);

	while ((*ACQ_IF_DUMP) == 1 && (*ACQ_MODE_V2) == 1)
	{
		// Use 1ms peak detector
		f_Peak_Detector(CH[c], Cpu_Z_IFFT);

		if (CH[c].detect == true)
		{
			(*ACTIVECHANNEL)++;
			Cpu_code_phase[c] = CH[c].codephase;

			if (start_mode == false)
			{
				if (CH[c].PRN < 32)
				{
					CH[c + 1].PRN = CH[c].PRN + 1;
					c++;
				}
				else
				{
					(*ACQ_MODE_V2) = 0;
				}
			}
			else
			{
				c++;

				if (c == (*ACTIVECHANNEL))
				{
					(*ACQ_MODE_V2) = 0;
				}
			}
		}
		else
		{
			if (CH[c].PRN < 32)
			{
				CH[c].PRN++;
			}
			else
			{
				(*ACQ_MODE_V2) = 0;
			}
		}

		(*ACQ_IF_DUMP) = 0;	// reset flag
	}

	cufftDestroy(p_cufft1);
	cufftDestroy(p_cufft2);
	cufftDestroy(p_cuifft);
}


__global__ void kernel_Correlation_First(short* PRN, short* code_index, short* carr_temp, int* s_c, short* ifdata, float* dump, int* cuda_CA_code_GPS_2046, int* trk_CH_2046)
{
	__shared__ int cache_dump_I_E[threadsPerBlock];		// ���� �޸� �Ҵ� -> �� �����尡 �ջ��� ����� �����ϴµ� ���
	__shared__ int cache_dump_I_P[threadsPerBlock];
	__shared__ int cache_dump_I_L[threadsPerBlock];
	__shared__ int cache_dump_Q_E[threadsPerBlock];
	__shared__ int cache_dump_Q_P[threadsPerBlock];
	__shared__ int cache_dump_Q_L[threadsPerBlock];

	int thread_num = threadIdx.x + (blockIdx.x * blockDim.x);
	int CH = blockIdx.y;

	short IF_DATA;
	bool Dump_Sync_Check = false;
	int thread_code_index = threadIdx.x + (blockIdx.x * blockDim.x) + ((int)(CorrSizeLocal * 2) * CH);

	int temp_I_E, temp_I_P, temp_I_L, temp_Q_E, temp_Q_P, temp_Q_L;
	int code_e, code_p, code_l;
	int carr_I, carr_Q;

	temp_I_E = temp_I_P = temp_I_L = temp_Q_E = temp_Q_P = temp_Q_L = carr_I = carr_Q = 0.0;
	code_e = code_p = code_l = 0;


	if (thread_num > trk_CH_2046[CH * 2])
		Dump_Sync_Check = true;


	IF_DATA = ifdata[thread_num + (int)(CorrSizeLocal * 2) * CH] >> 8;

	if (Dump_Sync_Check == false)
	{
		code_e = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] + 1 + s_c[CH]) + 2046) % 2046];
		code_p = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] + s_c[CH]) + 2046) % 2046];
		code_l = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] - 1 + s_c[CH]) + 2046) % 2046];
	}
	else
	{
		code_e = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] + 1) + 2046) % 2046];
		code_p = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index]) + 2046) % 2046];
		code_l = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] - 1) + 2046) % 2046];
	}

	switch (carr_temp[thread_code_index])
	{
	case 0:
		carr_I = 1;
		carr_Q = 2;
		break;
	case 1:
		carr_I = 2;
		carr_Q = 1;
		break;
	case 2:
		carr_I = 2;
		carr_Q = -1;
		break;
	case 3:
		carr_I = 1;
		carr_Q = -2;
		break;
	case 4:
		carr_I = -1;
		carr_Q = -2;
		break;
	case 5:
		carr_I = -2;
		carr_Q = -1;
		break;
	case 6:
		carr_I = -2;
		carr_Q = 1;
		break;
	default:
		carr_I = -1;
		carr_Q = 2;
	}



	__syncthreads();

	if (thread_num <= trk_CH_2046[CH * 2 + 1])
	{
		temp_I_E = IF_DATA * carr_I * code_e;
		temp_I_P = IF_DATA * carr_I * code_p;
		temp_I_L = IF_DATA * carr_I * code_l;
		temp_Q_E = IF_DATA * carr_Q * code_e;
		temp_Q_P = IF_DATA * carr_Q * code_p;
		temp_Q_L = IF_DATA * carr_Q * code_l;
	}
	else
	{
		temp_I_E = 0;
		temp_I_P = 0;
		temp_I_L = 0;
		temp_Q_E = 0;
		temp_Q_P = 0;
		temp_Q_L = 0;
	}

	__syncthreads();

	cache_dump_I_E[threadIdx.x] = 0;
	cache_dump_I_P[threadIdx.x] = 0;
	cache_dump_I_L[threadIdx.x] = 0;
	cache_dump_Q_E[threadIdx.x] = 0;
	cache_dump_Q_P[threadIdx.x] = 0;
	cache_dump_Q_L[threadIdx.x] = 0;

	__syncthreads();

	cache_dump_I_E[threadIdx.x] = temp_I_E;
	cache_dump_I_P[threadIdx.x] = temp_I_P;
	cache_dump_I_L[threadIdx.x] = temp_I_L;
	cache_dump_Q_E[threadIdx.x] = temp_Q_E;
	cache_dump_Q_P[threadIdx.x] = temp_Q_P;
	cache_dump_Q_L[threadIdx.x] = temp_Q_L;

	__syncthreads();


	int i = blockDim.x / 2;
	while (i != 0)
	{
		if (threadIdx.x < i)
		{
			cache_dump_I_E[threadIdx.x] += cache_dump_I_E[threadIdx.x + i];
			cache_dump_I_P[threadIdx.x] += cache_dump_I_P[threadIdx.x + i];
			cache_dump_I_L[threadIdx.x] += cache_dump_I_L[threadIdx.x + i];
			cache_dump_Q_E[threadIdx.x] += cache_dump_Q_E[threadIdx.x + i];
			cache_dump_Q_P[threadIdx.x] += cache_dump_Q_P[threadIdx.x + i];
			cache_dump_Q_L[threadIdx.x] += cache_dump_Q_L[threadIdx.x + i];

		}
		__syncthreads();


		i /= 2;
	}


	if (threadIdx.x == 0)
	{
		dump[(CH * gridDim.x * 6) + blockIdx.x] = cache_dump_I_E[0];
		dump[(CH * gridDim.x * 6) + (gridDim.x * 1) + blockIdx.x] = cache_dump_I_P[0];
		dump[(CH * gridDim.x * 6) + (gridDim.x * 2) + blockIdx.x] = cache_dump_I_L[0];
		dump[(CH * gridDim.x * 6) + (gridDim.x * 3) + blockIdx.x] = cache_dump_Q_E[0];
		dump[(CH * gridDim.x * 6) + (gridDim.x * 4) + blockIdx.x] = cache_dump_Q_P[0];
		dump[(CH * gridDim.x * 6) + (gridDim.x * 5) + blockIdx.x] = cache_dump_Q_L[0];

	}

}

__global__ void kernel_Correlation2(short* PRN, short* code_index, short* carr_temp, short* ifdata, float* dump, int* cuda_CA_code_GPS_2046, int* trk_CH_2046, int* trk_count)
{
	__shared__ int cache_dump_I_E[threadsPerBlock];		// ���� �޸� �Ҵ� -> �� �����尡 �ջ��� ����� �����ϴµ� ���
	__shared__ int cache_dump_I_P[threadsPerBlock];
	__shared__ int cache_dump_I_L[threadsPerBlock];
	__shared__ int cache_dump_Q_E[threadsPerBlock];
	__shared__ int cache_dump_Q_P[threadsPerBlock];
	__shared__ int cache_dump_Q_L[threadsPerBlock];
	
	int thread_num = threadIdx.x + (blockIdx.x * blockDim.x);
	int CH = blockIdx.y;
	int shift_thread_num = thread_num + trk_CH_2046[2 * CH] + 1;
	

	//if (thread_num > (trk_CH_2046[2 * CH + 1] - trk_CH_2046[2 * CH]))	return;

	short IF_DATA;
	//bool Dump_Sync_Check = false;
	int thread_code_index = thread_num + ((int)CorrSizeLocal * 2 * CH);
	//int thread_code_index = shift_thread_num + ((int)CorrSizeLocal * 2 * CH);
	
	int test_I_E, test_I_P, test_I_L, test_Q_E, test_Q_P, test_Q_L;

	int temp_I_E, temp_I_P, temp_I_L, temp_Q_E, temp_Q_P, temp_Q_L;
	int code_e, code_p, code_l;
	int carr_I, carr_Q;
	//short IF_DATA = 0;

	temp_I_E = temp_I_P = temp_I_L = temp_Q_E = temp_Q_P = temp_Q_L = carr_I = carr_Q = 0.0;
	code_e = code_p = code_l = 0;


	//if (thread_num > trk_CH_2046[CH * 2])
	//	Dump_Sync_Check == true;

	IF_DATA = ifdata[thread_num + (int)(CorrSizeLocal * 2) * CH] >> 8;
	//IF_DATA = ifdata[shift_thread_num] >> 8;


	//if (Dump_Sync_Check == false)
	//{
	//	code_e = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] + 1 + s_c[CH] * 2) + 2046) % 2046];
	//	code_p = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] + s_c[CH] * 2) + 2046) % 2046];
	//	code_l = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] - 1 + s_c[CH] * 2) + 2046) % 2046];
	//}
	
	
	code_e = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] + 1) + 2046) % 2046];
	code_p = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index]) + 2046) % 2046];
	code_l = cuda_CA_code_GPS_2046[((PRN[CH]) - 1) * 2046 + ((code_index[thread_code_index] - 1) + 2046) % 2046];
	

	switch (carr_temp[thread_code_index])
	{
	case 0:
		carr_I = 1;
		carr_Q = 2;
		break;
	case 1:
		carr_I = 2;
		carr_Q = 1;
		break;
	case 2:
		carr_I = 2;
		carr_Q = -1;
		break;
	case 3:
		carr_I = 1;
		carr_Q = -2;
		break;
	case 4:
		carr_I = -1;
		carr_Q = -2;
		break;
	case 5:
		carr_I = -2;
		carr_Q = -1;
		break;
	case 6:
		carr_I = -2;
		carr_Q = 1;
		break;
	default:
		carr_I = -1;
		carr_Q = 2;
	}

	__syncthreads();

	if (thread_num <= trk_CH_2046[CH * 2 + 1])
	//if (shift_thread_num <= trk_CH_2046[CH * 2 + 1])
	{
		temp_I_E = IF_DATA * carr_I * code_e;
		temp_I_P = IF_DATA * carr_I * code_p;
		temp_I_L = IF_DATA * carr_I * code_l;
		temp_Q_E = IF_DATA * carr_Q * code_e;
		temp_Q_P = IF_DATA * carr_Q * code_p;
		temp_Q_L = IF_DATA * carr_Q * code_l;
	}
	else
	{
		temp_I_E = 0;
		temp_I_P = 0;
		temp_I_L = 0;
		temp_Q_E = 0;
		temp_Q_P = 0;
		temp_Q_L = 0;
	}

	__syncthreads();

	cache_dump_I_E[threadIdx.x] = 0;
	cache_dump_I_P[threadIdx.x] = 0;
	cache_dump_I_L[threadIdx.x] = 0;
	cache_dump_Q_E[threadIdx.x] = 0;
	cache_dump_Q_P[threadIdx.x] = 0;
	cache_dump_Q_L[threadIdx.x] = 0;

	__syncthreads();

	
	cache_dump_I_E[threadIdx.x] = temp_I_E;
	cache_dump_I_P[threadIdx.x] = temp_I_P;
	cache_dump_I_L[threadIdx.x] = temp_I_L;
	cache_dump_Q_E[threadIdx.x] = temp_Q_E;
	cache_dump_Q_P[threadIdx.x] = temp_Q_P;
	cache_dump_Q_L[threadIdx.x] = temp_Q_L;

	__syncthreads();

	
	int i = blockDim.x / 2;
	while (i != 0)
	{
		if (threadIdx.x < i)
		{
			cache_dump_I_E[threadIdx.x] += cache_dump_I_E[threadIdx.x + i];
			cache_dump_I_P[threadIdx.x] += cache_dump_I_P[threadIdx.x + i];
			cache_dump_I_L[threadIdx.x] += cache_dump_I_L[threadIdx.x + i];
			cache_dump_Q_E[threadIdx.x] += cache_dump_Q_E[threadIdx.x + i];
			cache_dump_Q_P[threadIdx.x] += cache_dump_Q_P[threadIdx.x + i];
			cache_dump_Q_L[threadIdx.x] += cache_dump_Q_L[threadIdx.x + i];
		}
		__syncthreads();
		i /= 2;
	}

	if (threadIdx.x == 0)
	{
		dump[(blockIdx.y * gridDim.x * 6) + blockIdx.x] = cache_dump_I_E[0];
		dump[(blockIdx.y * gridDim.x * 6) + (gridDim.x * 1) + blockIdx.x] = cache_dump_I_P[0];
		dump[(blockIdx.y * gridDim.x * 6) + (gridDim.x * 2) + blockIdx.x] = cache_dump_I_L[0];
		dump[(blockIdx.y * gridDim.x * 6) + (gridDim.x * 3) + blockIdx.x] = cache_dump_Q_E[0];
		dump[(blockIdx.y * gridDim.x * 6) + (gridDim.x * 4) + blockIdx.x] = cache_dump_Q_P[0];
		dump[(blockIdx.y * gridDim.x * 6) + (gridDim.x * 5) + blockIdx.x] = cache_dump_Q_L[0];

	}

}

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
)
{

	if(TRK_first == 1)
	{	
		dim3 grids_trk(First_trk_blocksPerGrid, ACTIVECHANNEL);
		kernel_Correlation_First << <grids_trk, threadsPerBlock >> > (PRN, code_index, carr_temp, s_c, ifdata, dump, cuda_CA_code_GPS_2046, trk_CH_2046);

	}
	else
	{
		dim3 grids_trk(trk_blocksPerGrid, ACTIVECHANNEL);
		kernel_Correlation2 << <grids_trk, threadsPerBlock >> > (PRN, code_index, carr_temp, ifdata, dump, cuda_CA_code_GPS_2046, trk_CH_2046, trk_count);
	}

}

void f_constant_setup(int setup_f_s, int setup_f_if_GPS, int setup_f_ca_GPS, unsigned int sample_size_1ms)
{
	cudaMemcpyToSymbol(Gpu_const_f_s, &setup_f_s, sizeof(int));
	cudaMemcpyToSymbol(Gpu_const_f_if_GPS, &setup_f_if_GPS, sizeof(int));
	cudaMemcpyToSymbol(Gpu_const_f_ca_GPS, &setup_f_ca_GPS, sizeof(int));
	cudaMemcpyToSymbol(Gpu_const_sample_size_1ms, &sample_size_1ms, sizeof(unsigned long));

}