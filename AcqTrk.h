#pragma once
#include "stdafx.h"
#include "DSP.h"
#include "math.h"
#include <iostream>
//#include <Eigen/Dense>
#include "Global.h"
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
//#include <complex.h>

//#define _USE_MATH_DEFINES

int getMaxInt(int* n, int size);
int getMinInt(int* n, int size);
double getMaxdouble(double* n, int size);
double getMaxdouble_test(double* n, int size, int* i);
int getMaxIndex(double* Peak_Ratio, int size);
int getMaxIndex_PeakRatio(double* Peak_Ratio, int size);
int getMaxIndex(unsigned long long int* Peak_Ratio, int size);
double getMindouble(double* n, int size);

void f_Peak_Detector(Channel_struct& f_CH, double* Z_IFFT);
void f_Peak_Detector_5ms(Channel_struct& f_CH, double* Z_IFFT);
void f_loop_filter(Channel_struct& f_CH, int DCO_Carr_bit, int DCO_Code_bit, double DCO_Carr_RESOL, double DCO_Code_RESOL, int f_s);	// FLL, DLL, PLL ����
void f_Data_Read(int ifdata, int TRK_first, int ACTIVECHANNEL, short* TRK_CH_IF_DATA, short* TRK_IN_IF_DATA, int* N_trk, bool* Trk_ready, int* remaining_data);
void f_Code_Count(Channel_struct& f_CH, int TRK_first, int* first_check, int lnN, int* check_2046, int* remaining_data, short* Cur_IF_2bytes_sign, FILE* pFile_IFData, short* TRK_IN_IF_DATA, short* TRK_CH_IF_DATA, int k);
void f_1ms_Save(Channel_struct& f_CH);
void f_1ms_Read(Channel_struct& f_CH);
void f_BitCheck(Channel_struct& f_CH, double temp_phase_IQ);
void f_Tracking_GPS(int ifdata, int* TRK_first, int ACTIVECHANNEL, short* TRK_CH_IF_DATA, short* TRK_IN_IF_DATA, int* N_trk, bool* Trk_ready, int* remaining_data, short* Cur_IF_2bytes_sign, FILE* pFile_IFData, int* Gpu_CA_code_GPS_2046, int* CorrSizeLocal_count);

