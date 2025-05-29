#include "stdafx.h"
#include "Global.h"
#include "Defines.h"


double READTIME = 0.01;		
double CORRTIME = 0.001;

int f_s = Sampling_Frequency; // 25MHz
unsigned int ReadSizeLocal = int((double)f_s * READTIME);		
unsigned int CorrSizeLocal = int((double)f_s * CORRTIME);		

int f_if_GPS = IF_Frequency;
int f_ca_GPS = 1023000;	// CA_code_rate _GPS	// 1.023MHz

int CA_code_GPS[32][1023] = { 0, };			// GPS CA �ڵ� ������ ���� �ʱ�ȭ (32, 1023chip)
int CA_code_GPS_2046[32][2046] = { 0, };		// GPS CA �ڵ� ������ ���� �ʱ�ȭ (32, 1023chip) [32][2046]


double DCO_Code_RESOL = 0;
double DCO_Carr_RESOL = 0;

int g_DCO_Code_bit = 0;
int g_DCO_Carr_bit = 0;

unsigned int sample_size_1ms = f_s / 1000;
unsigned int sample_size_5ms = f_s * 5 / 1000;

int Ch_num_GPS = CHANNEL_MAX_NUM;
int	Ch_num = Ch_num_GPS;

int next_PRN_temp_GPS = 0;

//const int doppler_bin_size = 2 * (Doppler_Range * 1000) / Doppler_Interval + 1;

double FFT_acq_threshold = ACQ_Threshold;		//�⺻: 6.0


#ifdef ACTIVE_GPS_CHANNEL
int Z_threshold = TRK_Threshold;		// GPS (���� Bit�� Ȱ���Ͽ� ������ ��� (>>8) )
#endif

short range_dopp_freq_GPS[doppler_bin_size] = { 0, };			// ���÷� ���ļ� �˻� �� ���� : 45
double range_shift_code_GPS[2046] = { 0, };		// �ڵ� �˻� �� ���� : 2046



#ifdef ACTIVE_GPS_CHANNEL

int Ch_PRN_cold_GPS[ACTIVE_GPS_CHANNEL];


int Ch_PRN_warm_GPS[12] = { 28,3,17,22,1,11,19,0,0,0,0,0 };

int Ch_dopp_freq_warm_GPS[12] = { -150, 823, 1782, -1086, -2488, -2639, 2571, 2727, -2904, -2691, 0, 0 };	//�ʱ� ���÷� ���ļ� PRN �� ���� ���

#else
int Ch_PRN_cold_GPS[12] = { 0, };
int Ch_PRN_warm_GPS[12] = { 0, };
int Ch_dopp_freq_warm_GPS[12] = { 0, };
#endif

int Ch_dopp_index_warm_GPS[12] = { 0, };										//�ʱ� dopp index ���

bool start_mode = false;													// start_mode -> false : cold start
																			// start_mode -> true : warm start

Channel_struct CH[ACTIVE_GPS_CHANNEL] = { 0, };
Channel_struct NAV_CH[ACTIVE_GPS_CHANNEL] = { 0, };
Navigation_struct NAV = { 0, };				// �׹� ��� ����ü ����

iustruc ionoutc = { 0, };				/* Ionospheric model & UTC parameters. */




// ������
char FILE_name[100];
FILE* DataLog_DBG;
FILE* DataLog_Acq;
FILE* DataLog_NAVI;
FILE* TEST;


int SIGNAL_Time;
int  NAV_count = 0;

unsigned int sample_count = 0;

int ui_check = 0;

// CPU ó�� ����� ���� ����(ACQ)
cufftDoubleComplex* Cpu_replica_carr_vec = NULL;
cufftDoubleComplex* Cpu_replica_code_vec = NULL;
cufftDoubleComplex* Cpu_fft_out_carr_vec = NULL;
cufftDoubleComplex* Cpu_fft_out_code_vec = NULL;
cufftDoubleComplex* Cpu_freq_conj_carr_code = NULL;
cufftDoubleComplex* Cpu_dump_IFFT = NULL;
double* Cpu_Z_IFFT = NULL;
short* Cpu_ACQ_PRN = NULL;

// GPU ���� ó�� ����� ���� ����(ACQ)
short* Gpu_PRN = NULL;							// PRN ���� ���޿� GPU �޸�
short* Gpu_in_ifdata = NULL;					//if������ �Է¿�
short* Gpu_dopp_freq = NULL;
cufftDoubleComplex* Gpu_replica_carr_vec = NULL;
cufftDoubleComplex* Gpu_replica_code_vec = NULL;
cufftDoubleComplex* Gpu_fft_out_carr_vec = NULL;
cufftDoubleComplex* Gpu_fft_out_code_vec = NULL;
cufftDoubleComplex* Gpu_freq_conj_carr_code = NULL;
cufftDoubleComplex* Gpu_dump_IFFT = NULL;
double* Gpu_Z_IFFT = NULL;

// CPU ó�� ����� ���� ����(TRK)
float* Cpu_out_partial_dump = NULL;
short* Cpu_trk_PRN = NULL;
short* Cpu_code_index = NULL;
short* Cpu_carr_temp = NULL;
int* Cpu_code_phase = NULL;
int* Cpu_trk_CH_2046 = NULL;
int* Cpu_trk_corr_size = NULL;
int* Cpu_s_c = NULL;

int* Cpu_trk_count = NULL;

// GPU ���� ó�� ����� ���� ����(TRK)
short* Gpu_trk_in_ifdata = NULL;
short* Gpu_trk_ch_ifdata = NULL;
short* Gpu_trk_PRN = NULL;
float* Gpu_out_partial_dump = NULL;
short* Gpu_carr_temp = NULL;
short* Gpu_code_index = NULL;
int* Gpu_code_phase = NULL;
int* Gpu_trk_CH_2046 = NULL;
int* Gpu_trk_corr_size = NULL;
int* Gpu_s_c = NULL;

// GPU ������ �޸�
int* Gpu_trk_count = NULL;

int* Nav_DCO_Carr = NULL;
int* Nav_DCO_Code = NULL;
unsigned long* Nav_Epoch_Count_20m = NULL;
unsigned long* Nav_Frame_init_pos = NULL;
unsigned long* Nav_Epoch_Count = NULL;
unsigned long* Nav_Bit_init_pos = NULL;
unsigned long* Nav_Code_Count = NULL;
