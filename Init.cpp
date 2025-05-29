#pragma once
#include "stdafx.h"
#include "Init.h"
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "cublas_v2.h"
#include <cufft.h>

void f_GNSS_SDR_Init(Channel_struct *f_CH, Navigation_struct &NAV)	 // ���ű� �ʱ�ȭ
{
	int i = 0;

	int CA_GPS[1023] = { 0, };					// GPS CA �ڵ� ������ ���� �ʱ�ȭ (1, 1023chip)
	int k = 0;
	int kk = 0;

	// GPS Code ����
	for (k = 1; k <= 32; k++)				// 32���� CA �ڵ� ����
	{
		f_PRNcode(k, CA_GPS);
		for (kk = 0; kk<1023; kk++) 
		{
			CA_code_GPS[k - 1][kk] = CA_GPS[kk];
			CA_code_GPS_2046[k - 1][2 * kk] = CA_GPS[kk];
			CA_code_GPS_2046[k - 1][2 * kk + 1] = CA_GPS[kk];
		}
	}

	Cpu_replica_carr_vec = (cufftDoubleComplex*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(cufftDoubleComplex));
	Cpu_replica_code_vec = (cufftDoubleComplex*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(cufftDoubleComplex));
	Cpu_fft_out_carr_vec = (cufftDoubleComplex*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(cufftDoubleComplex));
	Cpu_fft_out_code_vec = (cufftDoubleComplex*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(cufftDoubleComplex));
	Cpu_freq_conj_carr_code = (cufftDoubleComplex*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(cufftDoubleComplex));
	Cpu_dump_IFFT = (cufftDoubleComplex*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(cufftDoubleComplex));
	Cpu_Z_IFFT = (double*)calloc(CorrSizeLocal * doppler_bin_size * 32, sizeof(double));
	Cpu_ACQ_PRN = (short*)calloc(32, sizeof(short));

	cudaMalloc((void**)&Gpu_PRN, sizeof(short) * 32);
	cudaMalloc((void**)&Gpu_in_ifdata, sizeof(short) * CorrSizeLocal * 5);
	cudaMalloc((void**)&Gpu_dopp_freq, sizeof(short) * doppler_bin_size);
	cudaMalloc((void**)&Gpu_replica_carr_vec, sizeof(cufftDoubleComplex) * CorrSizeLocal * doppler_bin_size * 32);
	cudaMalloc((void**)&Gpu_replica_code_vec, sizeof(cufftDoubleComplex) * CorrSizeLocal * doppler_bin_size * 32);
	cudaMalloc((void**)&Gpu_fft_out_carr_vec, sizeof(cufftDoubleComplex) * CorrSizeLocal * doppler_bin_size * 32);
	cudaMalloc((void**)&Gpu_fft_out_code_vec, sizeof(cufftDoubleComplex) * CorrSizeLocal * doppler_bin_size * 32);
	cudaMalloc((void**)&Gpu_freq_conj_carr_code, sizeof(cufftDoubleComplex) * CorrSizeLocal * doppler_bin_size * 32);
	cudaMalloc((void**)&Gpu_dump_IFFT, sizeof(cufftDoubleComplex) * CorrSizeLocal * doppler_bin_size * 32);
	cudaMalloc((void**)&Gpu_Z_IFFT, sizeof(double) * CorrSizeLocal * doppler_bin_size * 32);

	Cpu_trk_PRN = (short*)calloc(CHANNEL_MAX_NUM, sizeof(short));
	Cpu_out_partial_dump = (float*)calloc(blocksPerGrid * 6 * CHANNEL_MAX_NUM, sizeof(float));
	Cpu_carr_temp = (short*)calloc(CorrSizeLocal * CHANNEL_MAX_NUM * 2, sizeof(short));
	Cpu_code_index = (short*)calloc(CorrSizeLocal * CHANNEL_MAX_NUM * 2, sizeof(short));
	Cpu_code_phase = (int*)calloc(CHANNEL_MAX_NUM, sizeof(int));
	Cpu_trk_CH_2046 = (int*)calloc(CHANNEL_MAX_NUM * 2, sizeof(int));
	Cpu_trk_corr_size = (int*)calloc(1, sizeof(int));
	Cpu_s_c = (int*)calloc(CHANNEL_MAX_NUM, sizeof(double));

	Nav_DCO_Carr = (int*)calloc(CHANNEL_MAX_NUM, sizeof(int));
	Nav_DCO_Code = (int*)calloc(CHANNEL_MAX_NUM, sizeof(int));
	Nav_Epoch_Count_20m = (unsigned long*)calloc(CHANNEL_MAX_NUM, sizeof(unsigned long));
	Nav_Frame_init_pos = (unsigned long*)calloc(CHANNEL_MAX_NUM, sizeof(unsigned long));
	Nav_Epoch_Count = (unsigned long*)calloc(CHANNEL_MAX_NUM, sizeof(unsigned long));
	Nav_Bit_init_pos = (unsigned long*)calloc(CHANNEL_MAX_NUM, sizeof(unsigned long));
	Nav_Code_Count = (unsigned long*)calloc(CHANNEL_MAX_NUM, sizeof(unsigned long));

	Cpu_trk_count = (int*)calloc(1, sizeof(int));


	cudaMalloc((void**)&Gpu_trk_in_ifdata, sizeof(short) * CorrSizeLocal * 10);
	cudaMalloc((void**)&Gpu_trk_ch_ifdata, sizeof(short) * CorrSizeLocal * 2 * CHANNEL_MAX_NUM);
	cudaMalloc((void**)&Gpu_trk_PRN, sizeof(short) * CHANNEL_MAX_NUM);
	cudaMalloc((void**)&Gpu_out_partial_dump, sizeof(float) * blocksPerGrid * 6 * CHANNEL_MAX_NUM);
	cudaMalloc((void**)&Gpu_carr_temp, sizeof(short) * CorrSizeLocal * CHANNEL_MAX_NUM * 2);
	cudaMalloc((void**)&Gpu_code_index, sizeof(short) * CorrSizeLocal * CHANNEL_MAX_NUM * 2);
	cudaMalloc((void**)&Gpu_code_phase, sizeof(int) * CHANNEL_MAX_NUM);
	cudaMalloc((void**)&Gpu_trk_CH_2046, sizeof(int) * CHANNEL_MAX_NUM * 2);
	cudaMalloc((void**)&Gpu_trk_corr_size, sizeof(int));
	cudaMalloc((void**)&Gpu_s_c, sizeof(int) * CHANNEL_MAX_NUM);

	cudaMalloc((void**)&Gpu_trk_count, sizeof(int));



	for (k = 0; k < doppler_bin_size; k++)
	{
		range_dopp_freq_GPS[k] = (short)(pow((double)-1, (double)k) * 500 * (floor((double)((k + 1) / 2.0)))); // -11*10^3:500:11*10^3; // (Hz) ���÷� ���ļ� ����

	}



	// GPS Setting
	for (k = 0; k<2046; k++)
	{
		range_shift_code_GPS[k] = 0.5*k; // 0:1/2:1023-1/2; // (chip) �ڵ� õ�� ����
	}


	int DCO_Carr_bit = CARR_NC0_BIT;		//�⺻: 27
	g_DCO_Carr_bit = DCO_Carr_bit;
	//GPS Setting
	int DCO_Carr_INC_IF_GPS = (int)((f_if_GPS*pow((double)2, (double)DCO_Carr_bit)) / f_s);	
	DCO_Carr_RESOL = f_s / pow((double)2, (double)DCO_Carr_bit);			

	int DCO_Code_bit = CODE_NCO_BIT;		//�⺻: 26
	g_DCO_Code_bit = DCO_Code_bit;
	//GPS Setting
	int DCO_Code_INC_IF_GPS = (int)((f_ca_GPS*pow((double)2, (double)DCO_Code_bit)) / f_s);
	DCO_Code_RESOL = f_s / pow((double)2, (double)DCO_Code_bit);			

	// cold start , warm start ��� ����
	start_mode = false;			//false : cold start

	int Ch_PRN[ACTIVE_GPS_CHANNEL] = { 0, };

	for (int i = 0; i < ACTIVE_GPS_CHANNEL; i++)
	{
		Ch_PRN_cold_GPS[i] = i + 1;
	}

	double FLL_BW = FLL_Bandwidth;	//10
	double FLL_w0 = FLL_BW / 0.53;
	double FLL_coeff1 = FLL_w0 * FLL_w0;
	double FLL_coeff2 = M_SQRT2*FLL_w0;

	double PLL_BW = PLL_Bandwidth;	//20
	double PLL_w0 = PLL_BW / 0.7845;
	double PLL_coeff1 = PLL_w0*PLL_w0*PLL_w0;
	double PLL_coeff2 = 1.1*PLL_w0*PLL_w0;
	double PLL_coeff3 = 2.4*PLL_w0;

	double DLL_BW = DLL_Bandwidth; //0.25
	double DLL_w0 = DLL_BW / 0.25;
;
	int s_c_FFT_peak_index = 0;
	
	for (k = 0; k<Ch_num; k++)
	{
		
		if (k >= 0 && k<Ch_num_GPS)
		{
			
			CH[k].SYS = 0;	// GPS �Ҵ�	//0 ~ 11�� ä��
			CH[k].detect = false;
			CH[k].code_index = 0;
			CH[k].n = 0;

			if (start_mode == true) 
			{
				CH[k].PRN = Ch_PRN_warm_GPS[k];
				CH[k].d_f = Ch_dopp_freq_warm_GPS[k];
				CH[k].s_c = range_shift_code_GPS[0];
				CH[k].d_f_index = Ch_dopp_index_warm_GPS[k];
				CH[k].d_f_init = Ch_dopp_freq_warm_GPS[k];
				
			}
			else 
			{
				
				CH[k].PRN = Ch_PRN_cold_GPS[k];	//Ori
				
				//CH[k].PRN = 1; //MS
				CH[k].d_f = range_dopp_freq_GPS[0];
				CH[k].s_c = range_shift_code_GPS[0];
				CH[k].d_f_index = 0;
				CH[k].d_f_init = 0;
			}

			CH[k].d_f_code = CH[k].d_f/1540.0;

			CH[k].used_data = 0;

			CH[k].DCO_Carr_INC = DCO_Carr_INC_IF_GPS;		// ä�κ� GNSS �ݼ��� ���ļ� -> loop filter�� ���� ����
			CH[k].DCO_Code_INC = DCO_Code_INC_IF_GPS;		// ä�κ� GNSS �ڵ� ���ļ� -> loop filter�� ���� ����
			CH[k].DCO_Carr_INC_IF = DCO_Carr_INC_IF_GPS;		// ä�κ� GNSS IF ���ļ� �Ҵ�
			CH[k].DCO_Code_INC_IF = DCO_Code_INC_IF_GPS;		// ä�κ� GNSS Code ���ļ� �Ҵ�
			CH[k].FLL_coeff1 = FLL_coeff1;
			CH[k].FLL_coeff2 = FLL_coeff2;
			CH[k].PLL_coeff1 = PLL_coeff1;
			CH[k].PLL_coeff2 = PLL_coeff2;
			CH[k].PLL_coeff3 = PLL_coeff3;
			CH[k].DLL_w0 = DLL_w0;
			
		}

		CH[k].SampleSize_1ms = (unsigned long)((double)f_s/1000.0);

		CH[k].phase_code = 0.0;
		CH[k].phase_carr = 0.0;
		CH[k].codephase = 0;
		CH[k].dump_finish = 0;
		CH[k].store_finish = 0;

		CH[k].ACQ_mode = 1; // 1:FFT ��ȣ ȹ��		

		CH[k].Dump_Sync_Check = false;	// true = Dump Sync �Ϸ�
		CH[k].Dump_finish = false;

		if (CH[k].PRN == 0)
		{

			CH[k].ACQ_Check = false;

		}
		else
		{

			CH[k].ACQ_Check = true;
		}

		CH[k].TRC_Check = false;
		CH[k].DCO_Code_Check = true;												

		CH[k].Bit_Check = false;

		for (int ii = 0; ii <= 9; ii++)
		{
			CH[k].GPS_word[ii] = 0;
		}
		CH[k].FLL_Discri_inv_start = false;
		CH[k].FLL_Discri_Count = 0;
		CH[k].Epoch_Count_20m = 0; // 20msec ���� ī����
		CH[k].Epoch_Count = 0;	// 1msec ���� ī����
		CH[k].Code_Count = 0;	// 1chip ���� ī����
		CH[k].Code_Phase = 0;	// NCO�κ����� 1chip �� ����
		CH[k].Code_Meas = 0;	// Code ����ġ
		CH[k].Frame_init_pos = 0;
		CH[k].Frame_init_pos_Check = false;
		
	}
	
	NAV.GPS_rec_time = 0;															

	NAV.rec_TIC = 0;																
	NAV.rec_TIC_Count = 0;															
	NAV.TIC_Interval_s = 0.1;	//0.1													
	NAV.PLOT_Interval_s = 0.1;	//0.01												
	NAV.TIC_Interval_s_GPS = NAV.TIC_Interval_s;									
	NAV.rec_TIC_Count_Interval = (int)floor((double)(f_s*NAV.TIC_Interval_s));		
	NAV.rec_PLOT_Count_Interval = (int)floor((double)(f_s*NAV.PLOT_Interval_s));	

	NAV.GPS_L1_CA_NVS = 0;				// (NVS = Number of Visible Satellite)
	for (i=0 ; i<12 ; i++)
		NAV.GPS_L1_CA_SN[i] = 0;		// GPS ��ȿ ���� ��ȣ �Ҵ�

	// �׹����� �Ķ���� ���� �� �ʱ�ȭ
	NAV.Navi_Start_Check = false;
	NAV.Navi_Count = 0;						
	NAV.pos_ecef.x = 0;		NAV.pos_ecef.y = 0;		NAV.pos_ecef.z = 0;		// ���ű� ��ġ [x y z]
	NAV.pos_ecef_del.x = 0;	NAV.pos_ecef_del.y = 0;	NAV.pos_ecef_del.z = 0;	// ���ű� ��ġ ���� ��� del_pos x,y,z;
	NAV.pos_llh.lat = 0;	NAV.pos_llh.lon = 0;	NAV.pos_llh.hgt = 0;	// ���ű� ��ġ [Latitude Longitude Height]	
	NAV.vel_neu.x = 0;		NAV.vel_neu.y = 0;		NAV.vel_neu.z = 0;


	NAV.GPS_time_init_sync_Check = false;


	sample_count = 0;


}
