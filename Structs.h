#pragma once
#include "stdafx.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include <iostream>
#include "Defines.h"


#pragma once


#define INVALID    0                      /* Almanac and ephemeris status. */
#define VALID      1                      /* Almanac and ephemeris status. */

#define GPSV        37
#define GPLength    (long)1023
#pragma once

#define LogicHigh   (char)-1
#define LogicLow    (char)+1

#define GPSV        37
#define GPLength    (long)1023

#define SignExtend(x,y) (long)(((long)((x)<<(31-(y))))>>(31-(y)))

#define GLONASS_Data_Decoding(x,y,z) (unsigned __int64)((   ((unsigned __int64)(x)<<(62-(y)))&0x3FFFFFFFFFFFFFFF   )>>(62-(z)))

#define Uint64_NH_20(x,y) (unsigned __int64)( (((unsigned __int64)(x)<<((64-20)-(y)))) >>(64-20))
// x = 추출할 데이터, y = word 상에서 앞쪽 bit, z = 추출할 데이터 bit 수
#define Data_Emit_32bit(x,y,z) (unsigned long)( ((unsigned long)(x)<<(1+(y)))>>(32-(z)) )
// Parity Check 용
#define rotl(x,y) (((x)<<(y))|((x)>>(32-(y))))

// 위성 위치 계산 파라미터 정의
#define GravConstant 3.986005E14
#define SECONDS_IN_WEEK 604800.0
#define SECONDS_IN_HALF_WEEK 302400.0
#define SECONDS_IN_DAY 86400.0
#define SECONDS_IN_HOUR 3600.0
#define WGS84oe 7.2921151467E-5

// 위성 위치 계산 파라미터 정의(BeiDou)
#define CGCS2000oe 7.2921150E-5 // [rad/s]

// 항법을 위한 파라미터 정의
#define Speed_of_light 299792458
#define REARTH 6371006.84
#define R2D (180.0/M_PI)
//////////////////////////////Structure////////////////////////////////////


typedef struct
{
	double x;                                            /* x co-ordinate. */
	double y;                                            /* y co-ordinate. */
	double z;                                            /* z co-ordinate. */
} cartstruc;

// Longitude, latitude and height position state vector.
typedef struct
{
	double lat;                                   /* Latitude co-ordinate. */
	double lon;                                  /* Longitude co-ordinate. */
	double hgt;                                     /* Height co-ordinate. */
} llhstruc;

// Spherical vector.
typedef struct
{
	double el;                                   /* Elevation co-ordinate. */
	double az;                                     /* Azimuth co-ordinate. */
	double hgt;                                     /* Height co-ordinate. */
} sphstruc;

// Ionospheric and UTC data.
typedef struct
{
	short vflg;                                      /* 0 = no valid data. */
													 /* Ionospheric model parameters. */
	double a0, a1, a2, a3;                         /* AFCRL alpha parameters. */
	double b0, b1, b2, b3;                          /* AFCRL beta parameters. */
													/* UTC conversion parameters. */
	double A0, A1;                      /* Coeffs for determining UTC time. */
	unsigned long tot;     /* Reference time for A0 & A1, sec of GPS week. */
	short dtls;                           /* Cumulative past leap seconds. */
	unsigned wnt;                    /* Current UTC reference week number. */
	unsigned wnlsf;           /* Week number when dtlsf becomes effective. */
	short dn;           /* Day of week (1-7) when dtlsf becomes effective. */
	short dtlsf;                    /* Scheduled future leap second event. */
} iustruc;

struct Channel_struct {

	// 테스트 flag
	bool test_flag;

	unsigned long SampleSize_1ms;


	// Dump_Sync 테스트
	bool Dump_Sync_Check;	// true = Dump Sync 완료
							// false = Dump Sync 미완료
	bool Dump_finish;		// true = Dump finish 완료
							// false = Dump finish 미완료
							
	int SYS;	// 0 : GPS, 2 : BeiDou
	int PRN;
	int ACQ_mode;	// 0 : Serial Search, 1 : FFT Search
	unsigned short STATUS;	// #define CH_STATUS_NOTUSED		0x0000
							// #define CH_STATUS_USED			0x0001
							// #define CH_STATUS_ACQ			0x0002
							// #define CH_STATUS_TRK			0x0004
							// #define CH_STATUS_LOCK_BIT		0x0010
							// #define CH_STATUS_LOCK_HEADER	0x0020
							// #define CH_STATUS_LOCK_FRAME		0x0040
							// #define CH_STATUS_USED_EPH		0x0080
							// #define CH_STATUS_NAV			0x0100
	bool detect;
	double d_f;
	double s_c;
	double d_f_code;				// code doppler
	int DCO_Carr_INC;
	int DCO_Code_INC;
	int DCO_Carr_INC_IF;	// Carr Intermidiate frequency as GNSS
	int DCO_Code_INC_IF;	// Code as GNSS
	int d_f_index;
	int d_f_init;
	int s_c_index;
	int DCO_Carr;
	int DCO_Code;
	int code_index;
	bool Acq_Detect;
	bool ACQ_Check;
	bool TRC_Check;
	bool DCO_Code_Check;
	bool DCO_Code_Check_TRC;
	int n;


	// GPU 병렬 처리 연산용 위상 값 저장
	float phase_code;
	float phase_carr;

	int codephase;
	int dump_finish;
	int store_finish;

	int trk_2046_count;

	int onems_count;
	int onems_ready;
	int code_index_ready;

	int valid_data;	// 채널 버퍼에 남은 데이터
	int corr_size;	// 채널 별 corr 사이즈
	int used_data;	// 데이터 버퍼에서 사용한 데이터 개수

	double* Z_IFFT;
	int s_c_FFT_peak_index;

	int First_2046;
	int Second_2046;


	double carr_replica_I;
	double carr_replica_Q;
	int code_replica_E;
	int code_replica_P;
	int code_replica_L;

	int dump_count;
	int TRK_IF_count;

	int bf_carr;
	int bf_code;
	unsigned long bf_code_count;
	unsigned long bf_ephoch_count;
	unsigned long bf_ephoch_count_20;
	int bf_code_index;

	float temp_dump_I_E;
	float temp_dump_Q_E;
	float temp_dump_I_P;
	float temp_dump_Q_P;
	float temp_dump_I_L;
	float temp_dump_Q_L;
	float temp_dump_I_P_pre;
	float temp_dump_Q_P_pre;

	int Z;
	double DLL_Discri;
	double PLL_Discri;
	double FLL_Discri;
	int temp_dump_I_P_pre_PLL;
	int temp_dump_Q_P_pre_PLL;

	// Loop filter
	double FLL_BW;
	double FLL_coeff1;
	double FLL_coeff2;
	double FLL_A0[2];
	double FLL_A1[2];
	double FLL_temp;

	double DLL_BW;
	double DLL_w0;
	double DLL_A0[2];
	double DLL_A1[2];
	double DLL_A2[2];
	double DLL_temp;

	double PLL_BW;
	double PLL_coeff1;
	double PLL_coeff2;
	double PLL_coeff3;
	double PLL_n;
	double PLL_A0[2];
	double PLL_A1[2];
	double PLL_A2[2];
	double PLL_temp;

	double PLL_B0[3];
	double PLL_1to2;
	double PLL_B1[3];
	double PLL_2to3;

	double Filtered_FLL_Discri;
	double Filtered_DLL_Discri;
	double Filtered_PLL_Discri;

	int FLL_Discri_Count;		// 초기 위상차에 따른 주파수 변화 방향 유지를 위한 반전시점 지연_카운터
	bool FLL_Discri_inv_start;	// 초기 위상차에 따른 주파수 변화 방향 유지를 위한 반전시점 지연_플레그

	int Z_5;			// Z_5 : Tracking 유지 여부 판단용 파라미터
	int Z_5_index;		// Z_5 : 1msec X 5 횟수 Index


	// Bit Sync
	int Bit_Change_Count[20];
	int Bit_Change_Count_index_20msec;		//20msec를 위한 index
	int Bit_Change_Count_index_1000msec;	//1sec를 위한 index
	int Bit_Change_Location;
	bool Bit_Check;

	int Navi_Data[2000];	// 300(word) X 5(subframe) = 1500 [bit]
	int Navi_Data_Parity[2];// ((이전 word의 29,30bit) - Parity 목적) : 제거 예정
	int Navi_Data_30[30];	//							 Parity 목적 : 제거 예정
	int Navi_Data_300[300];	// 300Bit 단위로 데이터 int형태로 저장 ------------------------> 이후 GPS_word로 변환 (미 제거) 
	int Navi_Data_1500[1500];//							word 단위로 변환 : 제거 예정

	int Navi_Data_pre;		// 이전 시점(index)의 데이터 : 현재 시점의 데이터 산출용 
	int Navi_Data_now;		// 현재 시점(index)의 데이터
	int Navi_Data_index;

	// Frame Sync 
	bool Preamble_Check;	// true : Preamble_find
	bool Bit_Inversion;		// false : No_Inversion, true : Inversion
	int Subframe_ID_Check;	// Subframe ID = 1~5
	bool Zero_bits_Check;
	bool Parity_Check;
	int Parity_Fail_Count;	// Parity 실패 시 count
	bool Frame_Sync_Check;

	// GPS_word 데이터 관리용
	unsigned long Navi_Data_30bit[2];		// 초기 30bit 단위 데이터 저장용
	unsigned long GPS_word[10];			// word1~10:[0~9] // unit32
	unsigned long GPS_word_full[10][5]; // word1~10:[0~9], subframe1~5:[0~4] // unit32
	unsigned long GPS_word_index;		// word index = 0~29 = word1~word30

										// Data Decoding 
	bool Subframe_123_Check;		// Subframe_123_check : true (Subframe ID 1,2,3번 데이터 보유)
	bool Subframe_Input_Check[5];	// Subframe_Input_Check[0] : true (Subframe ID 1 데이터 보유) 
									// Subframe_Input_Check[1] : true (Subframe ID 2 데이터 보유)
									// Subframe_Input_Check[2] : true (Subframe ID 3 데이터 보유)
									// Subframe_Input_Check[3] : true (Subframe ID 4 데이터 보유)
									// Subframe_Input_Check[4] : true (Subframe ID 5 데이터 보유)

									// DataDecoding Parameter //
	unsigned long IODC;
	unsigned long Subframe_IODE[3];		// Subframe_IODE[0]	: Subframe ID 1의 IODE Number
										// Subframe_IODE[1]	: Subframe ID 2의 IODE Number
										// Subframe_IODE[2]	: Subframe ID 3의 IODE Number
	int IODE_pre;	// 이전 IODE number
	int IODE;		// 현재 IODE number

					// 임시저장용 변수
	unsigned long temp_navi_data;
	unsigned long temp_navi_data1;
	unsigned long temp_navi_data2;
	long temp_navi_data_sign;
	//--Subframe ID 1--//
	unsigned long WN;
	double S_A;
	double S_H;
	double T_GD;
	double t_oe;
	double a_f2;
	double a_f1;
	double a_f0;

	//--Subframe ID 2--//
	double C_rs;
	double Del_n;
	double M_0;
	double C_uc;
	double e_s;
	double C_us;
	double root_a_s;
	double t_oe2;

	//--Subframe ID 3--//
	double C_ic;
	double Ome_e;
	double C_is;
	double i_0;
	double C_rc;
	double ome;
	double Ome_dat;
	double idot;


	unsigned long Epoch_Count_2m;
	unsigned long Epoch_Count_20m;
	unsigned long Epoch_Count;
	unsigned long Code_Count;	
	double Code_Phase;	
	double Carrier_Phase;
	double Code_delay_time;
	double Code_Meas;	
	unsigned long Frame_init_pos; 
	unsigned long Bit_init_pos;	
	bool Frame_init_pos_Check;	
	double sverc;		

	cartstruc SV_Pos;

	cartstruc SV_Vel;		// 위성 속도 [x y z]
	cartstruc SV_Acc;		// 위성 가속도 [x y z]
	double SV_Clk_Offset;	// 위성 클럭 offset
	double SV_Clk_Drift;	// 위성 클럭 drift
	double SV_Sec;

	unsigned long Z_Count;


	double Pred_Range;		// 위성위치-사용자 위치간의 거리(매 LMS마다 업데이트됨)
	double Pred_delay_time; // 위성위치-사용자 위치간의 전파지연 시간(매 LMS 마다 업데이트됨)

	double Peak1;
	double Peak2;
	double Peak_Ratio;

	double Doppler;

	double EL;	
	double AZ;

	int M_count;
	double I_sum;
	double Q_sum;
	double IQS_sum;
	double PNM;
	int K_count;
	double avgPNM;
	double CN0;

	
};



struct Navigation_struct {

	double GPS_rec_time;								

	double rec_TIC;										
	int rec_TIC_Count;									
	double TIC_Interval_s;								
	double PLOT_Interval_s;								
	double TIC_Interval_s_GPS;							
	int rec_TIC_Count_Interval;							
	int rec_PLOT_Count_Interval;						
																				

	unsigned long GPS_L1_CA_NVS;	// (NVS = Number of Visible Satellite)
	unsigned long GPS_L1_CA_SN[12];

	bool Navi_Start_Check;
	int Navi_Count;		
	cartstruc pos_ecef;		// 수신기 위치 [x y z]
	cartstruc pos_ecef_del;	// 수신기 위치 추정 결과 del_pos x,y,z;
	llhstruc pos_llh;		// 수신기 위치 [Latitude Longitude Height]
	llhstruc Mean_pos_llh;
	
	cartstruc vel_ecef;
	cartstruc vel_neu;
	double Clk_Drift;
										
	bool GPS_time_init_sync_Check;			// false : GPS 초기 시각 동기 미수행
											// true : GPS 초기 시각 동기 완료



};
