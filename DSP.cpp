#include "stdafx.h"
#include "DSP.h"
#include "math.h"
#include <iostream>
#include <Eigen/Dense>	// Matrix 사용을 위한 Eigen 라이브러리
#include "Global.h"

using Eigen::MatrixXd;

#define _USE_MATH_DEFINES


double VectorLength(double x, double y, double z)
{
	return(sqrt(x*x + y*y + z*z));
}


void XYZToLatLonHgt(cartstruc *xyz, llhstruc *llh)
{
	double a;                                    /* WGS-84 semimajor axis. */
	double b;                            /* WGS-84 semiminor (polar) axis. */
	double e;                                /* WGS-84 first eccentricity. */
	double e2;                       /* WGS-84 first eccentricity squared. */
	double eps;                                      /* Convergence limit. */
	double dpi2;                                     /* Convergence limit. */
	double n,d,nph,rho,latold,hgtold;

	eps = 1.0E-13;
	dpi2 = 1.570796326794897E0;

	a = 6378137.0E0;
	b = 6356752.3142E0;
	e = 0.0818191908426E0;
	e2 = 0.00669437999013E0;

	/* If the position is on the Z axis, the computation must be handled as a special case to keep it from blowing up. */
	rho = VectorLength(xyz->x,xyz->y,0.0);
	if (rho <= eps)
	{
		llh->lat = dpi2;                 /* Come here if we are on the Z axis. */
		if(xyz->z<0.0)
			llh->lat = -llh->lat;
		llh->lon = 0.0E0;
		llh->hgt = fabs(xyz->z) - b;
		return;
	}

	/* Come here in the typical case.  Since latitude and spheroid height depend on one another, the solution must be done iteratively. */
	llh->lat = atan2(xyz->z,rho);
	llh->lon = atan2(xyz->y,xyz->x);
	llh->hgt = rho/cos(llh->lat);

	latold = llh->lat + 1.0E0;
	hgtold = llh->hgt + 1.0E0;

	while(fabs(llh->lat - latold)>=eps || fabs(llh->hgt-hgtold)>=0.01)
	{
		/* Require latitude to converge to about the precision of the machine, and height to the precision of about a centimeter. */
		latold = llh->lat;
		hgtold = llh->hgt;
		d = e*sin(latold);
		n = a/sqrt(1.0E0 - d*d);
		llh->hgt = rho/cos(latold) - n;
		nph = n + llh->hgt;
		d = 1.0E0 - e2*(n/nph);
		llh->lat = atan2(xyz->z,rho*d);
	}
}

// Variables for GPS C/A Code generation
char G1[GPLength], G2[GPLength];
long G2Delay[GPSV] = {
	5,   6,   7,   8,   17,  18,  139, 140, 141, 251,
	252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
	473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
	861, 862, 863, 950, 947, 948, 950
};
long Phase_Assign_G2[37][2] = {
	{1,3},	{1,4},	{1,5},	{1,6},	{1,8},	{1,9},	{1,10},	{1,11},
	{2,7},	{3,4},	{3,5},	{3,6},	{3,8},	{3,9},	{3,10},	{3,11},
	{4,5},	{4,6},	{4,8},	{4,9},	{4,10},	{4,11},	{5,6},	{5,8},
	{5,9},	{5,10},	{5,11},	{6,8},	{6,9},	{6,10},	{6,11},	{8,9},
	{8,10},	{8,11},	{9,10},	{9,11},	{10,11}
};



void f_PRNcode(long sv, int CA[1023])
{
	long i, j;
	char SR[10], G;

	// Generate GPS G1 code
	for( i=0; i<10; i++ )   SR[i] = LogicHigh;
	for( i=0; i<GPLength; i++ )
	{
		G1[i] = SR[9];
		G     = SR[2] * SR[9];  // G = 1 + X^3 + X^10

		// Shifting
		for( j=9; j>0; j-- )    SR[j] = SR[j-1];
		SR[0] = G;
	}

	// Generate GPS G2 code
	for( i=0; i<10; i++ )   SR[i] = LogicHigh;
	for( i=0; i<GPLength; i++ )
	{
		G2[i] = SR[9];          // G = 1 + X^2 + X^3 + X^6 + X^8 + X^9 + X^10
		G     = SR[1] * SR[2] * SR[5] * SR[7] * SR[8] * SR[9];

		// Shifting
		for( j=9; j>0; j-- )    SR[j] = SR[j-1];
		SR[0] = G;
	}

	// Form C/A code by multiplying G1 and G2
	for( i=0; i<GPLength; i++ )
	{
		// the delay index of G2 code
		j = (i+GPLength-G2Delay[sv-1]) % GPLength;

		CA[i] = (G1[i] * G2[j]);
	}
}

void f_Carr_DCO (int& DCO_Carr, int DCO_Carr_bit, int DCO_Carr_INC)
{
	DCO_Carr += DCO_Carr_INC;
	if (DCO_Carr > (int)((0x00000001<<DCO_Carr_bit) -1))
	{
		DCO_Carr -= (0x00000001<<DCO_Carr_bit);
	}
	
}

void f_Code_DCO (int& DCO_Code, int DCO_Code_bit, int DCO_Code_INC, bool& DCO_Code_Check, int PRN)
{
	DCO_Code += DCO_Code_INC;
	if (DCO_Code > (int)((0x00000001<<DCO_Code_bit) -1))
	{
		DCO_Code -= (0x00000001<<DCO_Code_bit);
		DCO_Code_Check = true;
	}
}


void f_GPS_Word0_Word1_check(unsigned long word0, unsigned long word1, bool &Frame_Sync_Check, Channel_struct &f_CH)
{

	unsigned long preamble, zero_bit;

	// zero bit 확인 및 반전 처리
	if(word1&0x00000001L)
	{
		word0 = ~word0;
		word1 = ~word1;
	}

	if(word0&0x40000000L)
		word0 ^= 0x3FFFFFC0L;
	if(word1&0x40000000L)
		word1 ^= 0x3FFFFFC0L;
	
	preamble = (Data_Emit_32bit(word0,1,8) & 0xFF);
	if(preamble == 0x8B){
		f_CH.Preamble_Check = true;
	}										

	if(f_CH.Preamble_Check == true){


		f_CH.Subframe_ID_Check = Data_Emit_32bit(word1,20,3);

		zero_bit = Data_Emit_32bit(word1,29,2) & 0x3;
		if(zero_bit == 0x00)
		{
			f_CH.Zero_bits_Check = true;
		}

		f_CH.Parity_Check = false;

		f_Parity_Check (word0, f_CH.Parity_Check);
		if(f_CH.Parity_Check == true)
		{
			f_Parity_Check (word1, f_CH.Parity_Check);
		}
												
	}

	if(f_CH.Parity_Check == true && f_CH.Zero_bits_Check == true &&
		f_CH.Subframe_ID_Check >= 1 && f_CH.Subframe_ID_Check <= 5)
	{
		Frame_Sync_Check = true;
	}
	else
	{
		Frame_Sync_Check = false;
	}
}


void f_Parity_Check (unsigned long word_now, bool &Parity_Check)
{
	unsigned long d1,d2,d3,d4,d5,d6,d7,t,parity;

	d1 = word_now & 0xFBFFBF00L;
	d2 = rotl(word_now,1) & 0x07FFBF01L;
	d3 = rotl(word_now,2) & 0xFC0F8100L;
	d4 = rotl(word_now,3) & 0xF81FFE02L;
	d5 = rotl(word_now,4) & 0xFC00000EL;
	d6 = rotl(word_now,5) & 0x07F00001L;
	d7 = rotl(word_now,6) & 0x00003000L;

	t = d1 ^ d2 ^ d3 ^ d4 ^ d5 ^ d6 ^ d7;

	parity = t ^ rotl(t,6) ^ rotl(t,12) ^ rotl(t,18) ^ rotl(t,24);
	parity = parity & 0x3F;

	if(parity == (word_now & 0x0000003F))
	{
		Parity_Check = true;
	}
	else
	{
		Parity_Check = false;
	}

}

void f_Data_Decoding (Channel_struct &f_CH)
{
	int word_index;
	for (word_index=0 ; word_index<=9 ; word_index++){						
		f_CH.GPS_word_full[word_index][f_CH.Subframe_ID_Check-1] = f_CH.GPS_word[word_index];
	}

	if (f_CH.Subframe_123_Check == false){
		switch (f_CH.Subframe_ID_Check){
		case 1:
			f_CH.Subframe_Input_Check[0] = true;
			break;
		case 2:
			f_CH.Subframe_Input_Check[1] = true;
			break;
		case 3:
			f_CH.Subframe_Input_Check[2] = true;
			break;
		default:
			break;
		}

		if(f_CH.Subframe_Input_Check[0] == true && f_CH.Subframe_Input_Check[1] == true && f_CH.Subframe_Input_Check[2] == true ){	// Subframe 1,2,3번 데이터가 모두 있는지 판단
			f_CH.Subframe_123_Check = true;
		}	
	}

	if (f_CH.Subframe_123_Check == true){
		f_CH.Subframe_IODE[0] = (f_CH.GPS_word_full[7][0]>>22)&0xFF;		// Subframe ID 1의 IODE (8번word : 1~8bit)
		f_CH.Subframe_IODE[1] = (f_CH.GPS_word_full[2][1]>>22)&0xFF;		// Subframe ID 2의 IODE (3번word : 1~8bit)
		f_CH.Subframe_IODE[2] = (f_CH.GPS_word_full[9][2]>>22)&0xFF;		// Subframe ID 3의 IODE (10번word : 1~8bit)

		if(f_CH.Subframe_IODE[0] == f_CH.Subframe_IODE[1] && f_CH.Subframe_IODE[0] == f_CH.Subframe_IODE[2]){

			int bit_index[2] = {0,0};
			
			// Subframe ID 1
			// WN 
			bit_index[0] = 1-1; bit_index[1] = 30-10;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.WN = (f_CH.temp_navi_data) * 1;							// 끝 bit
			
			// S_A
			bit_index[0] = 13-1; bit_index[1] = 30-4;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.S_A = double(f_CH.temp_navi_data);					// 끝 bit
			
			// S_H
			bit_index[0] = 17-1; bit_index[1] = 30-6;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.S_H = double(f_CH.temp_navi_data) * 1;					// 끝 bit
			
			// IODE
			bit_index[0] = 1-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[8-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.Subframe_IODE[0] = (f_CH.temp_navi_data);					// 끝 bit
			
			// T_GD
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[7-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,7);			//two's complement, with the sign bit occupying the MSB
			f_CH.T_GD = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-31);	// 끝 bit
			
			// t_oe
			bit_index[0] = 9-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[8-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.t_oe = double(f_CH.temp_navi_data) * pow((double)2,(double)4);	// 끝 bit
			
			// a_f2
			bit_index[0] = 1-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[9-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,7);			//two's complement, with the sign bit occupying the MSB
			f_CH.a_f2 = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-55);	// 끝 bit
			
			// a_f1
			bit_index[0] = 9-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[9-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.a_f1 = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-43);	// 끝 bit
			
			// a_f0
			bit_index[0] = 1-1; bit_index[1] = 30-22;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[10-1][1-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit  
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,21);			//two's complement, with the sign bit occupying the MSB
			f_CH.a_f0 = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-31);	// 끝 bit

			// Subframe ID 2
			// IODE
			bit_index[0] = 1-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.Subframe_IODE[1] = (f_CH.temp_navi_data);							// 끝 bit
			
			// C_rs
			bit_index[0] = 9-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.C_rs = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-5);	// 끝 bit
			
			// Del_n
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[4-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.Del_n = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-43) * M_PI;	// 끝 bit
			
			// M_0
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[4-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data1 = f_CH.temp_navi_data<<(32-8);				// 앞쪽 8bit로 이동
			
			// M_0
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[5-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data2 = f_CH.temp_navi_data;						// 뒤쪽 24bit로 이동
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data1 | f_CH.temp_navi_data2,31);			//two's complement, with the sign bit occupying the MSB
			f_CH.M_0 = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-31) * M_PI;
			
			// C_uc
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[6-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.C_uc = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-29);	// 끝 bit
			
			// e_s
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[6-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data1 = f_CH.temp_navi_data<<(32-8);				// 앞쪽 8bit로 이동
			
			// e_s
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[7-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data2 = f_CH.temp_navi_data;						// 뒤쪽 24bit로 이동
			f_CH.e_s = double(f_CH.temp_navi_data1 | f_CH.temp_navi_data2) * pow((double)2,(double)-33);

			// C_us
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[8-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.C_us = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-29);	// 끝 bit

			// root_a_s
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[8-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data1 = f_CH.temp_navi_data<<(32-8);				// 앞쪽 8bit로 이동
			// root_a_s
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[9-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data2 = f_CH.temp_navi_data;						// 뒤쪽 24bit로 이동
			f_CH.root_a_s = double(f_CH.temp_navi_data1 | f_CH.temp_navi_data2) * pow((double)2,(double)-19);
			
			// t_oe2
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[10-1][2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.t_oe2 = double(f_CH.temp_navi_data) * pow((double)2,(double)4);	// 끝 bit


			//// Subframe ID 3
			// C_ic
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.C_ic = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-29);	// 끝 bit

			// Ome_e
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[3-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data1 = f_CH.temp_navi_data<<(32-8);				// 앞쪽 8bit로 이동
			
			// Ome_e
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[4-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data2 = f_CH.temp_navi_data;						// 뒤쪽 24bit로 이동
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data1 | f_CH.temp_navi_data2,31);			//two's complement, with the sign bit occupying the MSB
			f_CH.Ome_e = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-31) * M_PI;

			// C_is
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[5-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.C_is = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-29);	// 끝 bit

			// i_0
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[5-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data1 = f_CH.temp_navi_data<<(32-8);				// 앞쪽 8bit로 이동
			
			// i_0
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[6-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data2 = f_CH.temp_navi_data;						// 뒤쪽 24bit로 이동
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data1 | f_CH.temp_navi_data2,31);			//two's complement, with the sign bit occupying the MSB
			f_CH.i_0 = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-31) * M_PI;

			// C_rc
			bit_index[0] = 1-1; bit_index[1] = 30-16;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[7-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,15);			//two's complement, with the sign bit occupying the MSB
			f_CH.C_rc = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-5);	// 끝 bit

			// ome
			bit_index[0] = 17-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[7-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data1 = f_CH.temp_navi_data<<(32-8);				// 앞쪽 8bit로 이동
			
			// ome
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[8-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data2 = f_CH.temp_navi_data;						// 뒤쪽 24bit로 이동
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data1 | f_CH.temp_navi_data2,31);			//two's complement, with the sign bit occupying the MSB
			f_CH.ome = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-31) * M_PI;

			// Ome_dat
			bit_index[0] = 1-1; bit_index[1] = 30-24;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[9-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,23);			//two's complement, with the sign bit occupying the MSB
			f_CH.Ome_dat = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-43) * M_PI;	// 끝 bit
			
			// IODE
			bit_index[0] = 1-1; bit_index[1] = 30-8;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[10-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.Subframe_IODE[2] = f_CH.temp_navi_data;	// 끝 bit

			// idot
			bit_index[0] = 9-1; bit_index[1] = 30-14;	
			f_CH.temp_navi_data = ((f_CH.GPS_word_full[10-1][3-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
			f_CH.temp_navi_data_sign = SignExtend(f_CH.temp_navi_data,13);			//two's complement, with the sign bit occupying the MSB
			f_CH.idot = double(f_CH.temp_navi_data_sign) * pow((double)2,(double)-43) * M_PI;	// 끝 bit
	
#if EPH_PLOT

			// Ephemeris 디버깅
			char FILE_name_GPS[100];
			FILE* DataLog_NAVI_GPS;
			sprintf_s(FILE_name_GPS, "DataLog\\DataLog_NAVI_GPS_Eph%d.txt", f_CH.PRN);
			fopen_s(&DataLog_NAVI_GPS, FILE_name_GPS, "a");
			fprintf(DataLog_NAVI_GPS, "\n-------------Subframe ID 1-------------");
			fprintf(DataLog_NAVI_GPS, "\nWN : %d", f_CH.WN);							// WN 
			fprintf(DataLog_NAVI_GPS, "\nS_A : %f", f_CH.S_A);							// S_A
			fprintf(DataLog_NAVI_GPS, "\nS_H : %f", f_CH.S_H);							// S_H
			fprintf(DataLog_NAVI_GPS, "\nIODE_Sub1 : %d", f_CH.Subframe_IODE[0]);		// IODE
			fprintf(DataLog_NAVI_GPS, "\nT_GD : %e", f_CH.T_GD);						// T_GD
			fprintf(DataLog_NAVI_GPS, "\nt_oe : %e", f_CH.t_oe);						// t_oe
			fprintf(DataLog_NAVI_GPS, "\na_f2 : %e", f_CH.a_f2);						// a_f2
			fprintf(DataLog_NAVI_GPS, "\na_f1 : %e", f_CH.a_f1);						// a_f1
			fprintf(DataLog_NAVI_GPS, "\na_f0 : %e", f_CH.a_f0);						// a_f0

			// Subframe ID 2
			fprintf(DataLog_NAVI_GPS, "\n-------------Subframe ID 2-------------");
			fprintf(DataLog_NAVI_GPS, "\nIODE_Sub2 : %d", f_CH.Subframe_IODE[1]);		// IODE
			fprintf(DataLog_NAVI_GPS, "\nC_rs : %e", f_CH.C_rs);						// C_rs
			fprintf(DataLog_NAVI_GPS, "\nDel_n : %e", f_CH.Del_n);						// Del_n
			fprintf(DataLog_NAVI_GPS, "\nM_0 : %e", f_CH.M_0);							// M_0
			fprintf(DataLog_NAVI_GPS, "\nC_uc : %e", f_CH.C_uc);						// C_uc
			fprintf(DataLog_NAVI_GPS, "\ne_s : %e", f_CH.e_s);							// e_s
			fprintf(DataLog_NAVI_GPS, "\nC_us : %e", f_CH.C_us);						// C_us
			fprintf(DataLog_NAVI_GPS, "\nroot_a_s : %e", f_CH.root_a_s);				// root_a_s
			fprintf(DataLog_NAVI_GPS, "\nt_oe2 : %e", f_CH.t_oe2);						// t_oe2

			//// Subframe ID 3
			fprintf(DataLog_NAVI_GPS, "\n-------------Subframe ID 3-------------");
			fprintf(DataLog_NAVI_GPS, "\nC_ic : %e", f_CH.C_ic);						// C_ic
			fprintf(DataLog_NAVI_GPS, "\nOme_e : %e", f_CH.Ome_e);						// Ome_e
			fprintf(DataLog_NAVI_GPS, "\nC_is : %e", f_CH.C_is);						// C_is
			fprintf(DataLog_NAVI_GPS, "\ni_0 : %e", f_CH.i_0);							// i_0
			fprintf(DataLog_NAVI_GPS, "\nC_rc : %e", f_CH.C_rc);						// C_rc
			fprintf(DataLog_NAVI_GPS, "\nome : %e", f_CH.ome);							// ome
			fprintf(DataLog_NAVI_GPS, "\nOme_dat : %e", f_CH.Ome_dat);					// Ome_dat
			fprintf(DataLog_NAVI_GPS, "\nIODE_Sub3 : %d", f_CH.Subframe_IODE[2]);		// IODE
			fprintf(DataLog_NAVI_GPS, "\nidot : %e", f_CH.idot);						// idot
			fprintf(DataLog_NAVI_GPS, "\n");
			fclose(DataLog_NAVI_GPS);
#endif	
			
		}else{
			
		}
		
		ionoutc.vflg = INVALID;

	}
}


void f_Sat_Pos (Channel_struct &f_CH, double GPS_rev_time)	// 위성 위치 계산 함수
{
	double to;                         /* Reference time of the data set. */
	double toc;                      /* Reference time of clock data set. */
	double deltan;              /* Mean motion delta from computed value. */
	double cuc;          /* Cosine harmonic correction to orbital radius. */
	double cus;            /* Sine harmonic corr to argument of latitude. */
	double cic;                   /* Cosine harmonic corr to inclination. */
	double cis;                     /* Sine harmonic corr to inclination. */
	double crc;          /* Cosine harmonic correction to orbital radius. */
	double crs;            /* Sine harmonic correction to orbital radius. */
	double ecc;                                          /* Eccentricity. */
	double sqrta;                       /* Square root of semimajor axis. */
	double m0;                                    /* Mean anomaly at TOE. */
	double om0;                                /* Right ascension at TOE. */
	double in0;                                    /* Inclination at TOE. */
	double olc;                            /* Argument of perigee at TOE. */
	double omd;                               /* Rate of right ascension. */
	double idot;                                  /* Rate of inclination. */
	double af0;               /* Satellite clock offset wrt GPS time (s). */
	double af1;              /* Satellite clock drift wrt GPS time (s/s). */
	double af2;        /* Satellite clock drift rate wrt GPS time (s/s2). */
	double tgd;           /* Group delay correction for the SV clock (s). */
	double relativistic; /* Relativistic correction for the SV clock (s). */

	double tk;                             /* The prediction interval (s). */
	double a;                                  /* The Semi-Major Axis (m). */
	double mk;                           /* The Mean Anomlay at tk (rads). */
	double mkdot;                                     /* Derivative of mk. */
	double ek;                      /* The Eccentric Anomlay at tk (rads). */
	double ekold;                 /* Used in the iterative solution of ek. */
	double ekdot;                                     /* Derivative of ek. */
	double cek,sek;                                /* cos(ek) and sin(ek). */
	double pk;                   /* The Argument of Latitude at tk (rads). */
	double pkdot;                                     /* Derivative of pk. */
	double c2pk,s2pk;                          /* cos(2*pk) and sin(2*pk). */
	double uk;         /* The corrected Argument of Latitude at tk (rads). */
	double ukdot;                                     /* Derivative of uk. */
	double cuk,suk;                                /* cos(uk) and sin(uk). */
	double ok;        /* The Longitude of the Ascending Node at tk (rads). */
	double okdot;                                     /* Derivative of ok. */
	double sok,cok;                                /* cos(ok) and sin(ok). */
	double ik;                      /* The Orbit Inclination at tk (rads). */
	double ikdot;                                     /* Derivative of ik. */
	double sik,cik;                                /* cos(ik) and sin(ik). */
	double rk;                            /* The Orbital Radius at tk (m). */
	double rkdot;                                     /* Derivative of rk. */
	double xpk,ypk;          /* Satellite vector in the orbital plane (m). */
	double xpkdot,ypkdot;                   /* Derivatives of xpk and ypk. */
	
	double mu_div_rk3, OneMinusecosE, sqrtOneMinuse2, tmp;

	// Ephemeris 데이터 로드

	to		= (double)f_CH.t_oe2; 
	toc		= (double)f_CH.t_oe; 
	deltan	= (double)f_CH.Del_n;
	cuc		= (double)f_CH.C_uc;
	cus		= (double)f_CH.C_us;
	cic		= (double)f_CH.C_ic;
	cis		= (double)f_CH.C_is;
	crc		= (double)f_CH.C_rc;
	crs		= (double)f_CH.C_rs;
	ecc		= (double)f_CH.e_s;		
	sqrta	= (double)f_CH.root_a_s;
	m0		= (double)f_CH.M_0;
	om0		= (double)f_CH.Ome_e;  
	in0		= (double)f_CH.i_0;
	olc		= (double)f_CH.ome;		
	omd		= (double)f_CH.Ome_dat;
	idot	= (double)f_CH.idot;
	af0		= (double)f_CH.a_f0;
	af1		= (double)f_CH.a_f1;
	af2		= (double)f_CH.a_f2;
	tgd		= (double)f_CH.T_GD;
	
	/* Compute the difference between the data reference time and the prediction time. */

	tk = GPS_rev_time - to;

	if(tk>SECONDS_IN_HALF_WEEK)
		tk -= SECONDS_IN_WEEK;
	else if(tk<-SECONDS_IN_HALF_WEEK)
		tk += SECONDS_IN_WEEK;
	
	/* Calculate the Semi-Major Axis (m). */
	a = sqrta*sqrta;

	/* Mean anomaly (rads) and its derivative (rads/s) at tk. */
	mkdot = sqrt(GravConstant/(a*a*a)) + deltan;
	mk = m0 + mkdot*tk;

	/*
	* Obtain the Eccentric Anomaly (rads) and its derivative (rads/sec) by solving Kepler's equation iteratively:
	*     M = E - e*sin(E)
	* Since the orbit is near circular, small eccentricity, the solution will converge rapidly.
	* The intial estimate of the Eccentric Anomaly is the Mean Anomaly.
	*/

	ek = mk;
	ekold = ek + 1.0;

	while(fabs(ek-ekold)>1.0E-14)
	{
		ekold = ek;
		OneMinusecosE = 1.0-ecc*cos(ekold);
		ek = ek + (mk-ekold+ecc*sin(ekold))/OneMinusecosE;
	}

	sek = sin(ek);
	cek = cos(ek);

	ekdot = mkdot/OneMinusecosE;

	/* Compute the relativistic correction term (seconds). */
	relativistic = -4.442807633E-10*ecc*sqrta*sek;

	/* Compute the argument of latitude (rads) and its derivative (rads/sec). */

	sqrtOneMinuse2 = sqrt(1.0 - ecc*ecc);

	pk = atan2(sqrtOneMinuse2*sek,cek-ecc) + olc;
	pkdot = sqrtOneMinuse2*ekdot/OneMinusecosE;

	s2pk = sin(2.0*pk);
	c2pk = cos(2.0*pk);

	/* Compute the corrected argument of latitude (rads) and its derivative (rads/sec). */
	uk = pk + cus*s2pk + cuc*c2pk;
	suk = sin(uk);
	cuk = cos(uk);
	ukdot = pkdot*(1.0 + 2.0*(cus*c2pk - cuc*s2pk));

	/* Compute the corrected radius (m) and its derivative (m/s). */
	rk = a*OneMinusecosE + crc*c2pk + crs*s2pk;
	rkdot = a*ecc*sek*ekdot + 2.0*pkdot*(crs*c2pk - crc*s2pk);

	/* Compute the corrected orbital inclination (rads) and its derivative (rads/s). */
	ik = in0 + idot*tk + cic*c2pk + cis*s2pk;
	sik = sin(ik);
	cik = cos(ik);
	ikdot = idot + 2.0*pkdot*(cis*c2pk - cic*s2pk);

	/* Compute the satellite's position vector in its orbital plane	and its derivative. */
	xpk = rk*cuk;
	ypk = rk*suk;
	xpkdot = rkdot*cuk - ypk*ukdot;
	ypkdot = rkdot*suk + xpk*ukdot;

	/* Compute the longitude of the ascending node (rads) and its derivative (rads/s). */

	okdot = omd - WGS84oe;
	ok = om0 + tk*okdot - WGS84oe*to;
	sok = sin(ok);
	cok = cos(ok);

	/* Compute the satellite's position in space. */

	f_CH.SV_Pos.x = xpk*cok - ypk*cik*sok; // X
	f_CH.SV_Pos.y = xpk*sok + ypk*cik*cok; // Y
	f_CH.SV_Pos.z = ypk*sik; // Z

	tmp = ypkdot*cik - ypk*sik*ikdot;

	f_CH.SV_Vel.x = -okdot*f_CH.SV_Pos.y + xpkdot*cok - tmp*sok;
	f_CH.SV_Vel.y = okdot*f_CH.SV_Pos.x + xpkdot*sok + tmp*cok;
	f_CH.SV_Vel.z = ypk*cik*ikdot + ypkdot*sik;

	mu_div_rk3 = - GravConstant/(rk*rk*rk);
	tmp = mu_div_rk3 + WGS84oe*WGS84oe;

	f_CH.SV_Acc.x = tmp*f_CH.SV_Pos.x + 2.0*f_CH.SV_Vel.y*WGS84oe;
	f_CH.SV_Acc.y = tmp*f_CH.SV_Pos.y - 2.0*f_CH.SV_Vel.x*WGS84oe;
	f_CH.SV_Acc.z = mu_div_rk3*f_CH.SV_Pos.z;

	/* Now calculate the current satellite clock offset, drift and drift rate with respect to GPS system time. */
	
	/* Compute the difference, in seconds, between toc and & sec. */
	 
	tk = GPS_rev_time - toc;

	if(tk > SECONDS_IN_HALF_WEEK)
		tk -= SECONDS_IN_WEEK;
	else if(tk < -SECONDS_IN_HALF_WEEK)
		tk += SECONDS_IN_WEEK;
	 

	f_CH.SV_Clk_Offset = af0 + tk*(af1 + tk*af2) + relativistic - tgd;  
	f_CH.SV_Clk_Drift = af1 + 2.0*tk*af2;  

}


void f_Pred_Range(Channel_struct &f_CH, cartstruc N_Pos)  // Pred_range 계산함수
{
	f_CH.Pred_Range = sqrt((f_CH.SV_Pos.x-N_Pos.x)*(f_CH.SV_Pos.x-N_Pos.x)
						+ (f_CH.SV_Pos.y-N_Pos.y)*(f_CH.SV_Pos.y-N_Pos.y)
						+ (f_CH.SV_Pos.z-N_Pos.z)*(f_CH.SV_Pos.z-N_Pos.z));
	f_CH.Pred_delay_time = f_CH.Pred_Range / Speed_of_light;
}

void f_Code_Meas(Channel_struct& f_CH, int DCO_Code_bit, double GPS_rev_time, cartstruc& N_Pos, Channel_struct CH)	// 코드 측정치 계산 함수
{

	double temp_20ms_cnt;
	double temp_1ms_cnt;
	double temp_1chip_cnt;

	double IONO_delay = 0;
	double TROP_delay = 0;
	llhstruc llh;
	
	f_CH.Carrier_Phase = (double)(f_CH.DCO_Carr) / (0x00000001 << 27);
	f_CH.Code_Phase = (double)(f_CH.DCO_Code) / (0x00000001 << DCO_Code_bit);
	
	temp_20ms_cnt = f_CH.Epoch_Count_20m + (f_CH.Frame_init_pos);
	temp_1ms_cnt = f_CH.Epoch_Count + (f_CH.Bit_init_pos);
	temp_1chip_cnt = ((double)f_CH.Code_Count / 2.0) + f_CH.s_c;			

	if (temp_1chip_cnt >= 1023) {
		temp_1chip_cnt = temp_1chip_cnt - 1023;
		temp_1ms_cnt = temp_1ms_cnt + 1;
	}
	if (temp_1ms_cnt >= 20) {
		temp_1ms_cnt = temp_1ms_cnt - 20;
		temp_20ms_cnt = temp_20ms_cnt + 1;
	}
	if (temp_20ms_cnt >= 50) {
		temp_20ms_cnt = temp_20ms_cnt - 50;
	}

	f_CH.Code_delay_time = temp_20ms_cnt * (0.02) + temp_1ms_cnt * (0.001) + (temp_1chip_cnt + f_CH.Code_Phase / 2) * (0.001 / 1023);
	
	double svatmos = 0;
	
	XYZToLatLonHgt(&N_Pos, &llh);
	Corrections(f_CH, GPS_rev_time, N_Pos, llh, IONO_delay, TROP_delay);

	// 지구 자전 오차 계산
	f_CH.sverc = WGS84oe / Speed_of_light * (f_CH.SV_Pos.y * (f_CH.SV_Pos.x - N_Pos.x)
		- f_CH.SV_Pos.x * (f_CH.SV_Pos.y - N_Pos.y));

	f_CH.Code_Meas = (fmod((double)(GPS_rev_time), (double)1) - f_CH.Code_delay_time);

	if (f_CH.Code_Meas > 0.5)
		f_CH.Code_Meas -= 1.0;
	else if (f_CH.Code_Meas < -0.5)
		f_CH.Code_Meas += 1.0;

	//f_CH.Code_Meas = (f_CH.Code_Meas + f_CH.SV_Clk_Offset) * Speed_of_light - f_CH.sverc - svatmos;

	f_CH.Code_Meas = f_CH.Code_Meas * (double)Speed_of_light;

#if Corr_SAT_CLK_EPH
	f_CH.Code_Meas += (f_CH.SV_Clk_Offset * (double)Speed_of_light);
#endif

#if Corr_EARTH_ROT
	f_CH.Code_Meas -= f_CH.sverc;
#endif

#if Corr_IONO
	f_CH.Code_Meas -= IONO_delay;
#endif

#if Corr_TROP
	f_CH.Code_Meas -= TROP_delay;
#endif

	

}


void f_Z_count_Upt (Channel_struct &f_CH)	// 현 시점의 frame의 Z-Count 추출
{

	int bit_index[2] = {0,0};				// data bit의 앞 및 뒤 shift bit 저장
			
	bit_index[0] = 1-1; bit_index[1] = 30-17;	
	f_CH.temp_navi_data = ((f_CH.GPS_word[2-1]<<bit_index[0]) & 0x3FFFFFFF)>>bit_index[1];	// 초기 bit 
	f_CH.Z_Count = unsigned long(f_CH.temp_navi_data) * 6;	// [Sec]

}

void f_Navi (Channel_struct *f_CH, cartstruc &N_Pos, cartstruc &N_del_Pos, double &GPS_rec_time, int Ch_num) // 항법 수행
{
	int k;				// 전체 채널 카운트
	int k_active = 0;	// 유효한 채널 카운트
	int Ch_num_active = 0;	// 유효한 채널 수

	// 항법을 위해 유효한 채널 개수 확인 -> 행렬 크기 결정
	for (k=0;k<Ch_num;k++){
		if((f_CH+k)->Subframe_123_Check == true){
			Ch_num_active++;
		}
	}

	// 항법 수행 행렬 선언
	MatrixXd H		= MatrixXd::Zero(Ch_num_active,4);	// H matrix
	MatrixXd del_x	= MatrixXd::Zero(4,1);				// del_x, del_y, del_z, delcB
	MatrixXd del_r	= MatrixXd::Zero(Ch_num_active,1);	// del_r(선형화 지점 거리 - 의사거리 측정치)

	k_active = 0;
	for (k=0;k<Ch_num;k++){
		if((f_CH+k)->Subframe_123_Check == true){
			del_r(k_active,0) = ((f_CH+k)->Pred_Range - (f_CH+k)->Code_Meas);
			k_active++;
		}
	}
	
	k_active = 0;
	for (k=0;k<Ch_num;k++){
		if((f_CH+k)->Subframe_123_Check == true){
			H(k_active,0) = ((f_CH+k)->SV_Pos.x - N_Pos.x)/(f_CH+k)->Pred_Range;
			H(k_active,1) = ((f_CH+k)->SV_Pos.y - N_Pos.y)/(f_CH+k)->Pred_Range;
			H(k_active,2) = ((f_CH+k)->SV_Pos.z - N_Pos.z)/(f_CH+k)->Pred_Range;
			H(k_active,3) = (1);
			k_active++;
		}
	}

	// 항법 수행
	del_x = (H.transpose() * H).inverse() * H.transpose() * del_r;
	
	// 위치, del 위치, 시각 업데이트
	N_del_Pos.x = (double)del_x(0,0);
	N_del_Pos.y = (double)del_x(1,0);
	N_del_Pos.z = (double)del_x(2,0);

	N_Pos.x += N_del_Pos.x;
	N_Pos.y += N_del_Pos.y;
	N_Pos.z += N_del_Pos.z;

	GPS_rec_time += (double)del_x(3,0) / Speed_of_light;	//del_x(3,0): 수신기 시계 오차
	
}


// Quantization_2bit
int Quantization(signed short ss, bool AGC_ON, unsigned int &cnt_Q, unsigned int &cnt_MSB_Q, double SAMPLING, double &Threshold2)
{
	int SIGN_MAG;
	double magnitude = 0;

	cnt_Q = cnt_Q + 1;

	// 2bit
	if(ss > Threshold2)
		SIGN_MAG = +3;
	else if(ss >= 0.0)
		SIGN_MAG = +1;
	else if(ss >= -Threshold2)
		SIGN_MAG = -1;
	else
		SIGN_MAG = -3;


	if (SIGN_MAG == 3 || SIGN_MAG == -3)	
		cnt_MSB_Q+=1;

	if(cnt_Q == (SAMPLING/1000) && AGC_ON == 1)
	{
		magnitude=(100.0*cnt_MSB_Q)/(SAMPLING/1000);

		if(magnitude<13.4)
			Threshold2-=(double)0.3;

		if(magnitude>13.4)
			Threshold2+=(double)0.3;

		cnt_Q = 0;
		cnt_MSB_Q = 0;
	}	

	return(SIGN_MAG);
}

void Corrections(Channel_struct &f_CH, double rev_time, cartstruc &N_Pos, llhstruc &N_Pos_llh, double &IONO, double &TROP)
{
	double elevation,azimuth;      /* The SV elevation and azimuth (rads). */
	double doppler;                                /* The SV Doppler (Hz). */
	double ionospheric_correction;      /* Correction to pseudo-range (m). */
	double tropospheric_correction;     /* Correction to pseudo-range (m). */
	double range;                        /* SV range WRT the observer. */
	double range_rate;              /* SV range_rate WRT the observer. */
	double t[3][3] = {0,};       /* XYZ to North-East-Down transformation matrix. */

	cartstruc xyz={0,};                /* SV range vector WRT the observer. */
	cartstruc neu={0,};      /* North-East-Up SV position WRT the observer. */


	/* Convert SV XYZ to NEU in observer's topocentric system. */
	xyz = SubVector(&f_CH.SV_Pos,&N_Pos);

	ECEFToTopoTransMat(&N_Pos,t);

	ECEFToTopo(&xyz,t,&neu);

	/* Calculate the SV elevation and azimuth. */
	elevation = atan2(neu.z,VectorLength(neu.x,neu.y,0.0));
	azimuth = atan2(neu.y,neu.x);

	if(azimuth<0.0)
		azimuth += (2.0*M_PI);
	
	f_CH.EL = elevation * R2D;
	f_CH.AZ = azimuth * R2D;


	/* Don't bother calculating ionospheric and tropospheric corrections when the satellite is below the horizon. */
	if(elevation < 0.0)
	{
		ionospheric_correction = 0.0;
		tropospheric_correction = 0.0;
	}
	else
	{
		double E;                     /* Observer longitude (semicircles). */
		double phi_u;                  /* Observer latitude (semicircles). */
		double lamda_u;               /* Observer longitude (semicircles). */
		double F;                                     /* Obliquity factor. */

		/* Translate the satellite elevation and the observer location into units of semi-circles. */
		E = elevation/M_PI;
		phi_u = N_Pos_llh.lat/M_PI;
		lamda_u = N_Pos_llh.lon/M_PI;

		F = 1.0 + 16.0*pow((0.53 - E),3.0);

		/* Use the nocturnal value when ionospheric corrections are	unavailable. */

		if(ionoutc.vflg!=VALID)
		{
			ionospheric_correction = F*5.0e-9*Speed_of_light;
		}
		else
		{
			double t;                              /* Observer local time. */
			double psi;              /* Earth central angle (semicircles). */
			double phi_i;       /* Sub-ionospheric latitude (semicircles). */
			double lamda_i;    /* Sub-ionospheric longitude (semicircles). */
			double phi_m;  /* Sub-iono geomagnetic latitude (semicircles). */
			double phi_m2;                                 /* phi_m*phi_m. */
			double phi_m3;                           /* phi_m*phi_m*phi_m. */
			double AMP;                   /* Vertical delay amplitude (s). */
			double PER;                        /* Period of the model (s). */
			double X;                        /* Phase of the model (rads). */
			double X2;                                             /* X*X. */
			double X4;                                           /* X2*X2. */

			ionoutc.a0 = 1.4901E-08;
			ionoutc.a1 = 2.2352E-08;
			ionoutc.a2 = -1.1921E-07;
			ionoutc.a3 = -1.1921E-07;

			ionoutc.b0 = 1.1264E+05;
			ionoutc.b1 = 1.3107E+05;
			ionoutc.b2 = -1.3107E+05;
			ionoutc.b3 = -1.9661E+05;


			psi = 0.0137/(E + 0.11) - 0.022;

			phi_i = phi_u + psi*cos(azimuth);


			if(phi_i>0.416)
				phi_i = 0.416;
			else if(phi_i<-0.416)
				phi_i = -0.416;

			lamda_i = lamda_u + psi*sin(azimuth)/cos(phi_i*M_PI);

			phi_m = phi_i + 0.064*cos((lamda_i-1.617)*M_PI);
			phi_m2 = phi_m*phi_m;
			phi_m3 = phi_m2*phi_m;

			AMP = ionoutc.a0 + ionoutc.a1*phi_m	+ ionoutc.a2*phi_m2 + ionoutc.a3*phi_m3;
			PER = ionoutc.b0 + ionoutc.b1*phi_m	+ ionoutc.b2*phi_m2 + ionoutc.b3*phi_m3;

			if(AMP<0.0)
				AMP = 0.0;

			if(PER<72000.0)
				PER = 72000.0;

			t = SECONDS_IN_DAY/2.0*lamda_i + rev_time;
			while(t>=SECONDS_IN_DAY)
				t -= SECONDS_IN_DAY;
			while(t<0)
				t += SECONDS_IN_DAY;

			X = 2.0*M_PI*(t - 50400.0)/PER;


			if(fabs(X)<1.57)
			{
				X2 = X*X;
				X4 = X2*X2;
				ionospheric_correction = F*(5.0e-9 + AMP*(1.0 - X2/2.0 + X4/24.0))*Speed_of_light;
			}
			else
				ionospheric_correction = F*5.0e-9*Speed_of_light;
		}

		tropospheric_correction = 2.4224*exp(-0.13346e-3*N_Pos_llh.hgt)/(0.026+sin(elevation));
		if(tropospheric_correction<0.0)
			tropospheric_correction = 0.0;
		else if (tropospheric_correction>200.0)
			tropospheric_correction = 200.0;
	}

	//svatmos = (float)ionospheric_correction + (float)tropospheric_correction;

	IONO = (float)ionospheric_correction;
	TROP = (float)tropospheric_correction;

}

void ECEFToTopo(cartstruc *ecef, double t[3][3], cartstruc *topo)
{
	topo->x = t[0][0]*ecef->x + t[0][1]*ecef->y + t[0][2]*ecef->z;
	topo->y = t[1][0]*ecef->x + t[1][1]*ecef->y + t[1][2]*ecef->z;
	topo->z = t[2][0]*ecef->x + t[2][1]*ecef->y + t[2][2]*ecef->z;
}

void ECEFToTopoTransMat(cartstruc *xyz, double t[3][3])
{
	sphstruc eah={0,};                      /* Observer spherical co-ordinates. */
	double sel, cel;                                  /* sin(el), cos(el). */
	double saz, caz;                                  /* sin(az), cos(az). */

	CartToSph(xyz,&eah);

	sel = sin(eah.el);
	cel = cos(eah.el);
	saz = sin(eah.az);
	caz = cos(eah.az);

	t[0][0] = -sel*caz;
	t[0][1] = -sel*saz;
	t[0][2] = cel;
	t[1][0] = -saz;
	t[1][1] = caz;
	t[1][2] = 0.0;
	t[2][0] = cel*caz;
	t[2][1] = cel*saz;
	t[2][2] = sel;
}


void CartToSph(cartstruc *xyz, sphstruc *eah)
{
	eah->el = atan2(xyz->z,VectorLength(xyz->x,xyz->y,0.0));
	eah->az = atan2(xyz->y,xyz->x);
	eah->hgt = VectorLength(xyz->x,xyz->y,xyz->z);
}


cartstruc SubVector(cartstruc *a, cartstruc *b)
{
	cartstruc c={0,};

	c.x = a->x - b->x;
	c.y = a->y - b->y;
	c.z = a->z - b->z;

	return(c);
}

void f_Navi_velocity(Channel_struct* f_CH, cartstruc& N_Pos, cartstruc& N_Vel, int Ch_num)
{
	int k;
	int k_active = 0;
	int Ch_num_active = 0;


	for (k = 0; k < Ch_num; k++)
	{
		if ((f_CH + k)->Subframe_123_Check == true)
		{
			Ch_num_active++;
		}
	}

	MatrixXd B = MatrixXd::Zero(Ch_num_active, 4);
	MatrixXd del_v = MatrixXd::Zero(4, 1);
	MatrixXd del_C = MatrixXd::Zero(Ch_num_active, 1);


	k_active = 0;
	for (k = 0; k < Ch_num; k++)
	{
		if ((f_CH + k)->Subframe_123_Check == true)
		{
			del_C(k_active, 0) = (((((f_CH + k)->SV_Pos.x - N_Pos.x) * (f_CH + k)->SV_Vel.x) + (((f_CH + k)->SV_Pos.y - N_Pos.y) * (f_CH + k)->SV_Vel.y) + (((f_CH + k)->SV_Pos.z - N_Pos.z) * (f_CH + k)->SV_Vel.z)) / (f_CH + k)->Pred_Range) + (L1_CYCLE * (f_CH + k)->Doppler);
			k_active++;
		}
	}


	k_active = 0;
	for (k = 0; k < Ch_num; k++)
	{
		if ((f_CH + k)->Subframe_123_Check == true)
		{
			B(k_active, 0) = ((f_CH + k)->SV_Pos.x - N_Pos.x) / (f_CH + k)->Pred_Range;
			B(k_active, 1) = ((f_CH + k)->SV_Pos.y - N_Pos.y) / (f_CH + k)->Pred_Range;
			B(k_active, 2) = ((f_CH + k)->SV_Pos.z - N_Pos.z) / (f_CH + k)->Pred_Range;
			B(k_active, 3) = (1);
			k_active++;
		}
	}

	del_v = (B.transpose() * B).inverse() * B.transpose() * del_C;


	N_Vel.x = (double)del_v(0, 0);
	N_Vel.y = (double)del_v(1, 0);
	N_Vel.z = (double)del_v(2, 0);

	NAV.Clk_Drift = (double)del_v(3, 0);

}

void XYZtoNED_Velocity(cartstruc* XYZ_vel, llhstruc* LLH_pos, cartstruc* NED_vel)
{

	double lat = LLH_pos->lat;
	double lon = LLH_pos->lon;
	double sin_lat = sin(lat);
	double cos_lat = cos(lat);
	double sin_lon = sin(lon);
	double cos_lon = cos(lon);

	double c_ne[3][3] =
	{
		{-sin_lat * cos_lon, -sin_lon, -cos_lat * cos_lon},
		{-sin_lat * sin_lon, cos_lon, -cos_lat * sin_lon},
		{cos_lat, 0, -sin_lat}
	};

	double c_en[3][3];

	for (int j = 0; j < 3; ++j)
	{
		for (int k = 0; k < 3; ++k)
		{
			c_en[j][k] = c_ne[k][j];
		}
	}

	double velned_temp[3];
	double velxyz[3];

	velxyz[0] = XYZ_vel->x;
	velxyz[1] = XYZ_vel->y;
	velxyz[2] = XYZ_vel->z;

	for (int j = 0; j < 3; ++j)
	{
		velned_temp[j] = 0;
		for (int k = 0; k < 3; ++k)
		{
			velned_temp[j] += c_en[j][k] * velxyz[k];
		}
	}

	double velned[3];

	for (int j = 0; j < 3; ++j)
	{
		velned[j] = velned_temp[j];
	}

	NED_vel->x = velned[0];
	NED_vel->y = velned[1];
	NED_vel->z = velned[2];

}

void CN0_Estimate(Channel_struct& f_CH, float I_P, float Q_P, int M_maximum, int K_maximum)
{

	if (f_CH.M_count < M_maximum)
	{
		f_CH.I_sum += abs(I_P);
		f_CH.Q_sum += abs(Q_P);
		f_CH.IQS_sum += (I_P * I_P + Q_P * Q_P);
		f_CH.M_count++;

		if (f_CH.M_count == M_maximum)
		{
			double PNK = f_CH.I_sum * f_CH.I_sum + f_CH.Q_sum * f_CH.Q_sum;
			double PWK = f_CH.IQS_sum;
			f_CH.PNM += (PNK / PWK);
			f_CH.K_count++;

			// Zero set
			f_CH.I_sum = 0;
			f_CH.Q_sum = 0;
			f_CH.IQS_sum = 0;
			f_CH.M_count = 0;

			if (f_CH.K_count == K_maximum)
			{
				f_CH.avgPNM = f_CH.PNM / (double)K_maximum;

				double cn0 = abs((1.0 / 0.001) * ((f_CH.avgPNM - 1.0) / ((double)M_maximum - f_CH.avgPNM)));

				f_CH.CN0 = 10.0f * (double)log10(cn0);
				f_CH.PNM = 0;
				f_CH.K_count = 0;
			}
		}

	}
}

void f_Navigation(int ACTIVECHANNEL)
{
	FILE* DataLog_NAVI_OUT = NULL;


	int k_active = 0;
	int k = 0;

	if (NAV.rec_TIC_Count == NAV.rec_TIC_Count_Interval)
	{
		NAV.rec_TIC_Count -= NAV.rec_TIC_Count_Interval;
		NAV.rec_TIC += NAV.TIC_Interval_s_GPS;
		NAV.GPS_rec_time += NAV.TIC_Interval_s_GPS;

		ui_check = 1;

		int jj = 0;
		// 가시 위성 수 카운트
		NAV.GPS_L1_CA_NVS = 0;

		for (jj = 0; jj < ACTIVECHANNEL; jj++)
		{
			// Subframe 1,2,3 체크를 통해 가시위성 입력
			if (NAV_CH[jj].Subframe_123_Check == true)
			{
				NAV.GPS_L1_CA_SN[NAV.GPS_L1_CA_NVS] = NAV_CH[jj].PRN;
				NAV.GPS_L1_CA_NVS++;
			}
		}

		if (NAV.GPS_L1_CA_NVS >= 4)	// GPS 위성 갯수가 4개 이상인 경우
		{

			if (NAV.GPS_time_init_sync_Check == false)
			{
 
				for (k = 0; k < ACTIVECHANNEL; k++)
				{
					if (NAV_CH[k].PRN == NAV.GPS_L1_CA_SN[0])
					{
						// Z-count update 
						f_Z_count_Upt(NAV_CH[k]);

						NAV.GPS_rec_time = (double)NAV_CH[k].Z_Count + 0.07; // 1980년 1월 5일 자정 기준의 (초)단위 값

						NAV.GPS_time_init_sync_Check = true;
						break;
					}
					else
					{
						NAV.GPS_time_init_sync_Check = false;
					}
				}

				// 항법 초기위치 설정
				if (NAV.GPS_time_init_sync_Check == true)
				{
					double range;
					double scale;
					NAV.pos_ecef.x = 0;
					NAV.pos_ecef.y = 0;
					NAV.pos_ecef.z = 0;
					k_active = 0;

					for (k = 0; k < ACTIVECHANNEL; k++)
					{
						if (NAV_CH[k].Subframe_123_Check == true)
						{
							f_Sat_Pos(NAV_CH[k], NAV.GPS_rec_time);	//현 시점 기준의 위성위치 계산
							NAV.pos_ecef.x += NAV_CH[k].SV_Pos.x;
							NAV.pos_ecef.y += NAV_CH[k].SV_Pos.y;
							NAV.pos_ecef.z += NAV_CH[k].SV_Pos.z;
							k_active += 1; // 유효 채널 카운트

						}
					}

					//평균 좌표 계산(수신기)
					NAV.pos_ecef.x = NAV.pos_ecef.x / k_active;
					NAV.pos_ecef.y = NAV.pos_ecef.y / k_active;
					NAV.pos_ecef.z = NAV.pos_ecef.z / k_active;
					range = VectorLength(NAV.pos_ecef.x, NAV.pos_ecef.y, NAV.pos_ecef.z);
					scale = REARTH / range;

					NAV.pos_ecef.x = NAV.pos_ecef.x * scale;
					NAV.pos_ecef.y = NAV.pos_ecef.y * scale;
					NAV.pos_ecef.z = NAV.pos_ecef.z * scale;

				}
			}
			if (NAV.GPS_time_init_sync_Check == true)
			{
				for (NAV.Navi_Count = 0; NAV.Navi_Count < 5; NAV.Navi_Count++)
				{
					for (k = 0; k < ACTIVECHANNEL; k++)
					{

						if (NAV_CH[k].Subframe_123_Check == true)
						{
							f_Sat_Pos(NAV_CH[k], (NAV.GPS_rec_time - (double)0.07));
							f_Pred_Range(NAV_CH[k], NAV.pos_ecef);

							f_Sat_Pos(NAV_CH[k], NAV.GPS_rec_time - NAV_CH[k].Pred_delay_time);	
							f_Pred_Range(NAV_CH[k], NAV.pos_ecef);

							f_Code_Meas(NAV_CH[k], g_DCO_Code_bit, NAV.GPS_rec_time, NAV.pos_ecef, NAV_CH[0]);

							NAV.Navi_Start_Check = true;

						}

					}

					if (NAV.Navi_Start_Check == true)
					{
						f_Navi(NAV_CH, NAV.pos_ecef, NAV.pos_ecef_del, NAV.GPS_rec_time, Ch_num);
						XYZToLatLonHgt(&NAV.pos_ecef, &NAV.pos_llh);

						f_Navi_velocity(NAV_CH, NAV.pos_ecef, NAV.vel_ecef, Ch_num);
						XYZtoNED_Velocity(&NAV.vel_ecef, &NAV.pos_llh, &NAV.vel_neu);

						// 항법 결과 출력
						if (NAV.Navi_Count == 4)
						{
							NAV_count++;
							NAV.Mean_pos_llh.lat = NAV.Mean_pos_llh.lat + NAV.pos_llh.lat;
							NAV.Mean_pos_llh.lon = NAV.Mean_pos_llh.lon + NAV.pos_llh.lon;
							NAV.Mean_pos_llh.hgt = NAV.Mean_pos_llh.hgt + NAV.pos_llh.hgt;

#if NAVI_PLOT		
							// 항법 결과 디버깅
							fopen_s(&DataLog_NAVI_OUT, "DataLog\\DataLog_Navi_Out.txt", "a");
							fprintf(DataLog_NAVI_OUT,
								"%f,%d,"
								"%f,%f,%f,"
								"%f,%f,%f,"
								"%f,%f,%f,"
								"%f,%f\n",
								NAV.rec_TIC, NAV.Navi_Count,
								NAV.pos_ecef.x, NAV.pos_ecef.y, NAV.pos_ecef.z,
								NAV.pos_llh.lat * (180.0 / M_PI), NAV.pos_llh.lon * (180.0 / M_PI), NAV.pos_llh.hgt,
								0, 0, 0,
								0, 0);
							//fprintf(DataLog_NAVI_OUT,
							//	"%f,%f,%d,"
							//	"%f,%f,%f,"
							//	"%f,%f,%f,"
							//	"%f,%f,%f,"
							//	"%f,%f\n",
							//	NAV.GPS_rec_time, NAV.rec_TIC, NAV.Navi_Count,
							//	NAV.pos_ecef.x, NAV.pos_ecef.y, NAV.pos_ecef.z,
							//	NAV.pos_llh.lat* (180.0 / M_PI), NAV.pos_llh.lon* (180.0 / M_PI), NAV.pos_llh.hgt,
							//	0, 0, 0,
							//	0, 0);
							fclose(DataLog_NAVI_OUT);

							// 항법 좌표 디버깅
							sprintf_s(FILE_name, "DataLog\\DataLog_Navi_LLH_Out.txt");
							fopen_s(&DataLog_NAVI_OUT, FILE_name, "a");
							fprintf(DataLog_NAVI_OUT, "%f	%f	%f	%f\n ", NAV.rec_TIC, NAV.pos_llh.lat * (180.0 / M_PI), NAV.pos_llh.lon * (180.0 / M_PI), NAV.pos_llh.hgt);
							fclose(DataLog_NAVI_OUT);

#endif

#if VEL_PLOT
							// 속도 결과 디버깅
							fopen_s(&DataLog_NAVI_OUT, "DataLog\\DataLog_Vel_Out.txt", "a");
							fprintf(DataLog_NAVI_OUT, "%f	%f	%f	%f\n", NAV.rec_TIC, NAV.vel_neu.x, NAV.vel_neu.y, NAV.vel_neu.z);
							fclose(DataLog_NAVI_OUT);
#endif		

#if ELAZ_PLOT
							// Elevation & Azimuth 디버깅
							if (NAV_count == 10)
							{
								for (int b = 0; b < ACTIVECHANNEL; b++)
								{
									sprintf_s(FILE_name, "DataLog\\DataLog_ELAZ_Out.txt");
									fopen_s(&DataLog_DBG, FILE_name, "a");
									fprintf(DataLog_DBG, "%d	%f	%f\n", NAV_CH[b].PRN, NAV_CH[b].EL, NAV_CH[b].AZ);
									fclose(DataLog_DBG);
								}
							}
#endif

						}
					}
				}
			}
		}
		else if (NAV.GPS_L1_CA_NVS <= 3)
		{
			NAV.GPS_time_init_sync_Check = false;
		}

	}
}