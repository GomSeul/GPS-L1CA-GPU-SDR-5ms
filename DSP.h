#pragma once
//#include "stdafx.h"
//#include "math.h"
//#include <iostream>
//#include <Eigen/Dense>
//#include <fftw3.h>
//#include "Global.h"
//
//#define _USE_MATH_DEFINES
//
//#pragma once

#include "Structs.h"


double VectorLength(double x, double y, double z);
void XYZToLatLonHgt(cartstruc *xyz, llhstruc *llh);

void f_PRNcode(long sv, int CA[1023]);

void f_Carr_DCO(int& DCO_Carr, int DCO_Carr_bit, int DCO_Carr_INC);
void f_Code_DCO(int& DCO_Code, int DCO_Code_bit, int DCO_Code_INC, bool& DCO_Code_Check, int PRN);


void f_Parity_Check(unsigned long word_now, bool &Parity_Check);
void f_GPS_Word0_Word1_check(unsigned long Word0, unsigned long Word1, bool &Frame_Sync_Check, Channel_struct &f_CH);

void f_Data_Decoding(Channel_struct &f_CH);

void f_Sat_Pos(Channel_struct &f_CH, double GPS_rec_time);

void f_Code_Meas(Channel_struct& f_CH, int DCO_Code_bit, double GPS_rev_time, cartstruc& N_Pos, Channel_struct CH);

void f_Z_count_Upt(Channel_struct &f_CH);
void f_Pred_Range(Channel_struct &f_CH, cartstruc N_Pos);	

void f_Navi(Channel_struct *f_CH, cartstruc &N_Pos, cartstruc &N_del_Pos, double &GPS_rec_time, int Ch_num);				
																																
int Quantization(signed short ss, bool AGC_ON, unsigned int &cnt_Q, unsigned int &cnt_MSB_Q, double SAMPLING, double &Threshold2);

void Corrections(Channel_struct& f_CH, double rev_time, cartstruc& N_Pos, llhstruc& N_Pos_llh, double& IONO, double& TROP);
void ECEFToTopo(cartstruc *ecef, double t[3][3], cartstruc *topo);
void ECEFToTopoTransMat(cartstruc *xyz, double t[3][3]);
void CartToSph(cartstruc *xyz, sphstruc *eah);
cartstruc SubVector(cartstruc *a, cartstruc *b);
void f_Navi_velocity(Channel_struct* f_CH, cartstruc& N_Pos, cartstruc& N_Vel, int Ch_num);
void XYZtoNED_Velocity(cartstruc* XYZ_vel, llhstruc* LLH_pos, cartstruc* NED_vel);

void CN0_Estimate(Channel_struct& f_CH, float I_P, float Q_P, int M_maximum, int K_maximum);

void f_Navigation(int ACTIVECHANNEL);
