#pragma once
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include "Defines.h"
#include "Structs.h"
#include <cufft.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"


extern double READTIME;		
extern double CORRTIME;		
extern int f_s;
extern unsigned int ReadSizeLocal;		
extern unsigned int CorrSizeLocal;		
extern int f_if_GPS;
extern int f_if_BDS;
extern int f_ca_GPS;
extern int f_ca_BDS;

extern int CA_code_GPS[32][1023];			
extern int CA_code_GPS_2046[32][2046];	


extern double DCO_Code_RESOL;
extern double DCO_Carr_RESOL;

// �������� code, carrier dco bit �߰�	20171224
extern int g_DCO_Code_bit;
extern int g_DCO_Carr_bit;

extern unsigned int sample_size_1ms;
extern unsigned int sample_size_5ms;

extern int Ch_num_GPS;
extern int Ch_num;

extern int next_PRN_temp_GPS;

//extern const int doppler_bin_size;
constexpr int doppler_bin_size = 2 * (Doppler_Range * 1000) / Doppler_Interval + 1;
extern double FFT_acq_threshold;

extern int Z_threshold;	

extern short range_dopp_freq_GPS[doppler_bin_size];
extern double range_shift_code_GPS[2046];	

extern int Ch_PRN_cold_GPS[ACTIVE_GPS_CHANNEL];
extern int Ch_PRN_warm_GPS[12];				//�ʱ� ä�� PRN ���� ��� �ʱ�ȭ warm
extern int Ch_dopp_freq_warm_GPS[12];		//�ʱ� ���÷� ���ļ� PRN �� ���� ���
extern int Ch_dopp_index_warm_GPS[12];		//�ʱ� dopp index ���

extern bool start_mode;										// start_mode -> false : cold start
															// start_mode -> true : warm start

extern Channel_struct CH[ACTIVE_GPS_CHANNEL];
extern Channel_struct NAV_CH[ACTIVE_GPS_CHANNEL];
extern Navigation_struct NAV;				// �׹� ��� ����ü ����

extern iustruc ionoutc;							/* Ionospheric model & UTC parameters. */

extern char FILE_name[100];
extern FILE *DataLog_DBG;
extern FILE *DataLog_Acq;
extern FILE *DataLog_NAVI;
extern FILE *TEST;

extern int SIGNAL_Time;
extern int NAV_count;

extern unsigned int sample_count;

extern int ui_check;

// CPU ó�� ����� ���� ����(ACQ)
extern cufftDoubleComplex* Cpu_replica_carr_vec;
extern cufftDoubleComplex* Cpu_replica_code_vec;
extern cufftDoubleComplex* Cpu_fft_out_carr_vec;
extern cufftDoubleComplex* Cpu_fft_out_code_vec;
extern cufftDoubleComplex* Cpu_freq_conj_carr_code;
extern cufftDoubleComplex* Cpu_dump_IFFT;
extern double* Cpu_Z_IFFT;
extern short* Cpu_ACQ_PRN;

// GPU ���� ó�� ����� ���� ����(ACQ)
extern short* Gpu_PRN;							// PRN ���� ���޿� GPU �޸�
extern short* Gpu_in_ifdata;					//if������ �Է¿�
extern short* Gpu_dopp_freq;
extern cufftDoubleComplex* Gpu_replica_carr_vec;
extern cufftDoubleComplex* Gpu_replica_code_vec;
extern cufftDoubleComplex* Gpu_fft_out_carr_vec;
extern cufftDoubleComplex* Gpu_fft_out_code_vec;
extern cufftDoubleComplex* Gpu_freq_conj_carr_code;
extern cufftDoubleComplex* Gpu_dump_IFFT;
extern double* Gpu_Z_IFFT;

// CPU ó�� ����� ���� ����(TRK)
extern float* Cpu_out_partial_dump;
extern short* Cpu_trk_PRN;
extern short* Cpu_code_index;
extern short* Cpu_carr_temp;
extern int* Cpu_code_phase;
extern int* Cpu_trk_CH_2046;
extern int* Cpu_trk_corr_size;
extern int* Cpu_s_c;

extern int* Cpu_trk_count;

// GPU ���� ó�� ����� ���� ����(TRK)
extern short* Gpu_trk_in_ifdata;
extern short* Gpu_trk_ch_ifdata;
extern short* Gpu_trk_PRN;
extern float* Gpu_out_partial_dump;
extern short* Gpu_carr_temp;
extern short* Gpu_code_index;
extern int* Gpu_code_phase;
extern int* Gpu_trk_CH_2046;
extern int* Gpu_trk_corr_size;
extern int* Gpu_s_c;

// GPU ������ �޸�
extern int* Gpu_trk_count;


extern int* Nav_DCO_Carr;
extern int* Nav_DCO_Code;
extern unsigned long* Nav_Epoch_Count_20m;
extern unsigned long* Nav_Frame_init_pos;
extern unsigned long* Nav_Epoch_Count;
extern unsigned long* Nav_Bit_init_pos;
extern unsigned long* Nav_Code_Count;


#endif