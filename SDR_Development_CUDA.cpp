// SDR_Development_CUDA_op.cpp : �ܼ� ���� ���α׷��� ���� �������� �����մϴ�.
//

#include "stdafx.h"
#include "CUDA_op.cuh"
#define _USE_MATH_DEFINES
#include "math.h"
#include "DSP.h"
#include "AcqTrk.h"
#ifdef _WIN32
    #include <Windows.h>
#endif
#include <iostream>
#include <Eigen/Dense>    // Matrix operations using Eigen library
#include "Global.h"
#include "Init.h"
#include "UI.h"
#include "stdio.h"
#include "Defines.h"

// For cross-platform large file support
#ifdef _WIN32
    #include <io.h>
    #define fseeko _fseeki64
    #define ftello _ftelli64
#else
    #define _FILE_OFFSET_BITS 64
#endif


using namespace std;

int main()
{
	// Console â ����
	// Console window setup
#ifdef _WIN32
	system("mode con:cols=180 lines=140");
	system("mode con:cols=180 lines=30");
	system("title GNSS SDR");
#else
	printf("GNSS SDR\n");
#endif

	FILE* DataLog_DBG = NULL;
	FILE* DataLog_NAVI_OUT;
	FILE* DataLog_CARR = NULL;

	// ä�� �ʱ�ȭ �� �Ķ���� �ʱ�ȭ �Լ�
	f_GNSS_SDR_Init(CH, NAV);	 // ���ű� �ʱ�ȭ

	// short* Cur_IF_2bytes_sign = NULL; // Not needed - reading directly into ifdata
	bool ACQ_operation = true;
	int N_acq = 0;
	int N_trk = 0;
	bool Trk_ready = false;

	FILE* pFile_IFData;

	int remaining_data = 0;
	int CorrSizeLocal_count = 0;

	int err;  // Using int instead of errno_t for cross-platform compatibility
	// Cross-platform file opening
	pFile_IFData = fopen("D:/GPS L1CA signal/GPS_L1CA_SIGNAL_CPU_IGS_100s_I_AMP_2_trajectory.bin", "rb");
	err = (pFile_IFData == NULL) ? 1 : 0;

	if (err != 0)
	{
		printf("file open error!\n");
#ifdef _WIN32
		system("pause");
#else
		printf("Press Enter to continue...\n");
		getchar();
#endif
		return 0;
	}

	
	// Cross-platform file seeking
	fseeko(pFile_IFData, 0, SEEK_END);
	unsigned long long SizeFull = ftello(pFile_IFData);
	SIGNAL_Time = (((SizeFull * 8) / ADC) / f_s);

	fseeko(pFile_IFData, 0, SEEK_SET);

	// Cur_IF_2bytes_sign = (short*)calloc(ReadSizeLocal, sizeof(short)); // Not needed

	short* ACQ_IF_DATA = NULL;
	ACQ_IF_DATA = (short*)calloc(CorrSizeLocal * 5, sizeof(short));

	short* TRK_IF_DATA = NULL;
	TRK_IF_DATA = (short*)calloc(CorrSizeLocal * CHANNEL_MAX_NUM, sizeof(short));

	short* TRK_IN_IF_DATA = NULL;
	TRK_IN_IF_DATA = (short*)calloc(CorrSizeLocal * 10, sizeof(short));

	short* TRK_CH_IF_DATA = NULL;	
	TRK_CH_IF_DATA = (short*)calloc((CorrSizeLocal * 2) * CHANNEL_MAX_NUM, sizeof(short));

	// GPU ��� �޸� ����
	f_constant_setup(f_s, f_if_GPS, f_ca_GPS, CorrSizeLocal);

	// GPU ���� ó�� ����� ���� ����
	size_t size_GPS_CA_CODE = sizeof(int) * 32 * 1023;
	int* Gpu_CA_code_GPS_1023 = NULL;
	cudaMalloc((void**)&Gpu_CA_code_GPS_1023, size_GPS_CA_CODE);
	cudaMemcpy(Gpu_CA_code_GPS_1023, CA_code_GPS, size_GPS_CA_CODE, cudaMemcpyHostToDevice);

	size_t size_GPS_CA_CODE_2046 = sizeof(int) * 32 * 2046;
	int* Gpu_CA_code_GPS_2046 = NULL;
	cudaMalloc((void**)&Gpu_CA_code_GPS_2046, size_GPS_CA_CODE_2046);
	cudaMemcpy(Gpu_CA_code_GPS_2046, CA_code_GPS_2046, size_GPS_CA_CODE_2046, cudaMemcpyHostToDevice);

	short ifdata = 0;								
	int ACTIVECHANNEL = 0;
	int ACQ_MODE_V2 = 1;
	int ACQ_IF_DUMP = 0;
	int TRK_first = 1;
	int init_time = 1;

	//1ms ���� ó�� ���� ����
	while (!feof(pFile_IFData))
	{
		fread(&ifdata, 1, sizeof(short), pFile_IFData);


		if (ACQ_MODE_V2 == 1)
		{
			ACQ_IF_DATA[N_acq] = ifdata;

			if (N_acq == CorrSizeLocal * 5)
			{
				ACQ_IF_DUMP = 1;	//5msec data read
				N_acq = 0; // Reset counter for next acquisition
			}

			N_acq++;

			if (ACQ_IF_DUMP == 1)
			{
				for (int count = 0; count < 32; count++)
				{
					Cpu_ACQ_PRN[count] = count + 1;
				}

				cudaMemcpy(Gpu_in_ifdata, ACQ_IF_DATA, sizeof(short) * CorrSizeLocal * 5, cudaMemcpyHostToDevice);
				cudaMemcpy(Gpu_dopp_freq, range_dopp_freq_GPS, sizeof(short) * doppler_bin_size, cudaMemcpyHostToDevice);
				cudaMemcpy(Gpu_PRN, Cpu_ACQ_PRN, sizeof(short) * 32, cudaMemcpyHostToDevice);

				f_Acquisition_5ms	// GPU Acquisition �Լ� ȣ��
				(
					Gpu_PRN,
					Gpu_in_ifdata,
					Gpu_CA_code_GPS_1023,
					Gpu_dopp_freq,
					Gpu_replica_carr_vec,
					Gpu_replica_code_vec,
					Gpu_fft_out_carr_vec,
					Gpu_fft_out_code_vec,
					Gpu_freq_conj_carr_code,
					Gpu_dump_IFFT,
					Gpu_Z_IFFT,
					Cpu_Z_IFFT,
					&ACQ_IF_DUMP,
					&ACQ_MODE_V2,
					&ACTIVECHANNEL,
					Cpu_code_phase
				);

				for (int i = 0; i < ACTIVECHANNEL; i++)
				{
					CH[i].onems_count = CorrSizeLocal;
				}

			}
			
		}


		// ��ȣ���� ����
		if (ACQ_MODE_V2 == 0) //ACQ �Ϸ�
		{
			f_Tracking_GPS
			(
				ifdata, 
				&TRK_first, 
				ACTIVECHANNEL, 
				TRK_CH_IF_DATA, 
				TRK_IN_IF_DATA, 
				&N_trk, 
				&Trk_ready, 
				&remaining_data, 
				&Cur_IF_2bytes_sign, 
				pFile_IFData, 
				Gpu_CA_code_GPS_2046, 
				&CorrSizeLocal_count
			);
		}
		
		// �׹� ����
		f_Navigation(ACTIVECHANNEL);
	
		// �ܼ�â ����
		f_Console_Plot_GPS(NAV_CH, NAV);

		if (feof(pFile_IFData) != 0)	// ���� ���� ������ ���
		{
			printf("file end -> processing end!\n");

	#ifdef _WIN32
		system("pause");
#else
		printf("Press Enter to continue...\n");
		getchar();
#endif
			return 0;

		}
	} // <----------------------- while (!feof(pFile_IFData))

	// GPU ���� �޸𸮸� �����Ѵ�.
	cudaFree(Gpu_PRN);
	cudaFree(Gpu_in_ifdata);
	cudaFree(Gpu_dopp_freq);
	cudaFree(Gpu_replica_carr_vec);
	cudaFree(Gpu_replica_code_vec);
	cudaFree(Gpu_CA_code_GPS_1023);

	cudaFree(Gpu_fft_out_carr_vec);
	cudaFree(Gpu_fft_out_code_vec);
	cudaFree(Gpu_freq_conj_carr_code);
	cudaFree(Gpu_dump_IFFT);
	cudaFree(Gpu_Z_IFFT);

	cudaFree(Gpu_trk_in_ifdata);
	cudaFree(Gpu_trk_PRN);
	cudaFree(Gpu_out_partial_dump);
	cudaFree(Gpu_carr_temp);
	cudaFree(Gpu_code_index);
	cudaFree(Gpu_code_phase);

	// CPU ���� �޸𸮸� �����Ѵ�.
	free(ACQ_IF_DATA);
	free(TRK_IF_DATA);
	free(TRK_IN_IF_DATA);
	free(Cpu_replica_carr_vec);
	free(Cpu_replica_code_vec);
	free(Cpu_fft_out_carr_vec);
	free(Cpu_fft_out_code_vec);
	free(Cpu_freq_conj_carr_code);
	free(Cpu_dump_IFFT);
	free(Cpu_Z_IFFT);

	free(Cpu_out_partial_dump);
	free(Cpu_carr_temp);
	free(Cpu_code_index);
	free(Cpu_code_phase);
	free(Cpu_ACQ_PRN);

	// ���� �ݱ�
	fclose(pFile_IFData);

	return 0;
}
