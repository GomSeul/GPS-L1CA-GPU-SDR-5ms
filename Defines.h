#pragma once

const int N = 25000;
const int N_5ms = 125000;  // 5ms sample size for acquisition
const int threadsPerBlock = 512;//32768 16384 8192 4096 2048 1024//512;		// threadPerBlock�� ������� ���� ��Ȯ�� ������, ���� ó�� �ּ� ��� ũ�� ����(���� �ӵ� ����)
const int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;	
const int blocksPerGrid_5ms = (N_5ms + threadsPerBlock - 1) / threadsPerBlock;  // Blocks for 5ms acquisition
const int First_trk_blocksPerGrid = (2 * N + threadsPerBlock - 1) / threadsPerBlock;
const int trk_blocksPerGrid = ((N + 100) + threadsPerBlock - 1) / threadsPerBlock;

// SDR ���� ����
#define Corr_EARTH_ROT		1
#define Corr_SAT_CLK_EPH	1
#define Corr_IONO			1
#define	Corr_TROP			1

// ����� ���� �Ķ����
#define Sampling_Frequency		25e6
#define IF_Frequency			5.42e6

#define Doppler_Range			11
#define Doppler_Interval		100

#define CARR_NC0_BIT			30
#define CODE_NCO_BIT			30	

#define FLL_Bandwidth			10			
#define PLL_Bandwidth			5			// 20
#define DLL_Bandwidth			0.25

#define ACQ_Threshold			5.5		
#define TRK_Threshold			1500

#define N_Addi_ACQ				3			// Origianl + Additional search Ƚ��

#define FLL_assisted_PLL		0			// �⺻ : FLL_ONLY




// MATLAB ������
#define ACQ_PLOT_1				0
#define ACQ_PLOT_2				0
#define TRK_PLOT				0
#define	EPH_PLOT				0
#define NAVI_PLOT				1
#define VEL_PLOT				0
#define ELAZ_PLOT				0



// CONSTANT VALUES
#define M_max					20
#define K_max					50

#define L1						1575.42e6
#define L1_CYCLE				0.190293672798365

#define D2R						(M_PI / 180.0)
#define N_Doppler_bin			(((Doppler_Range * 1000) / Doppler_Interval) * 2) + 1

#define ADC						16
#define CHANNEL_MAX_NUM			32
#define ACTIVE_GPS_CHANNEL		CHANNEL_MAX_NUM










