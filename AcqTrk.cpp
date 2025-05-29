#include "stdafx.h"
#include "AcqTrk.h"
#include "DSP.h"
#include "math.h"
#include <iostream>
#include "CUDA_op.cuh"

//#include <Eigen/Dense>	// Matrix ����� ���� Eigen ���̺귯��

//---------------Matrix ����� ���� �ڷ��� ����----------
//using Eigen::MatrixXd;


//#define _USE_MATH_DEFINES


float getMaxfloat(float* n, int size) {
	float max = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] > max) max = n[i];

	return max;
}

double getMaxdouble_test(double* n, int size, int* i)
{

	double max = n[0];
	
	for (int b = 1; b < size; b++)
		if (n[b] > max) 
		{
			max = n[b];
			*i = b; 
		}


	return max;
}

double getMaxdouble(double* n, int size)
{
	double max = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] > max) max = n[i];

	return max;


}

int getMaxIndex(double* Peak_Ratio, int size)
{
	unsigned long long int max = Peak_Ratio[0];

	int maxIndex = 0;

	for (int i = 1; i < size; i++)
	{
		if (Peak_Ratio[i] > max)
		{
			max = Peak_Ratio[i];
			maxIndex = i;
		}
	}
	return maxIndex;
}

int getMaxInt(int* n, int size) {
	int max = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] > max) max = n[i];

	return max;
}

int getMinInt(int* n, int size) {
	int min = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] < min) min = n[i];

	return min;
}

int getMaxIndex_PeakRatio(double* Peak_Ratio, int size)
{
	double max = Peak_Ratio[0];

	int maxIndex = 0;

	for (int i = 1; i < size; i++)
	{
		if (Peak_Ratio[i] > max)
		{
			max = Peak_Ratio[i];
			maxIndex = i;
		}
	}
	return maxIndex;
}

int getMaxIndex(unsigned long long int* Peak_Ratio, int size)
{
	unsigned long long int max = Peak_Ratio[0];

	int maxIndex = 0;

	for (int i = 1; i < size; i++)
	{
		if (Peak_Ratio[i] > max)
		{
			max = Peak_Ratio[i];
			maxIndex = i;
		}
	}
	return maxIndex;
}

double getMindouble(double* n, int size) {
	double min = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] < min) min = n[i];

	return min;
}

void f_Peak_Detector(Channel_struct& f_CH, double* Z_IFFT)
{

	int target_PRN = 10;
	int target_d_f_index = 11;

	double Peak1[45] = { 0 };
	double Peak2[45] = { 0 };
	int Peak_index[45] = { 0 };
	double Peak_Ratio[45] = { 0 };
	int s_c_FFT_peak_index[45] = { 0 };
	double* Z_IFFT_for_Peak2 = NULL;
	Z_IFFT_for_Peak2 = (double*)malloc(sizeof(double) * (int)(f_s / 1000.0));
	double* Z_IFFT_DUMP;
	Z_IFFT_DUMP = (double*)malloc(sizeof(double) * (int)(f_s / 1000.0));
	int detect_finish = 1;
	int MAX_INDEX = 0;

	while (f_CH.d_f_index < 45)
	{

		for (unsigned long jj = 0; jj < (sample_size_1ms); jj++)
		{
			Z_IFFT_DUMP[jj] = Z_IFFT[(CorrSizeLocal * 45 * (f_CH.PRN - 1)) + (CorrSizeLocal * f_CH.d_f_index) + jj];
		}

#if ACQ_PLOT_1

		sprintf_s(FILE_name, "DataLog\\Datalog_Acq_Full_FFT_GPS%d.txt", f_CH.PRN);
		fopen_s(&DataLog_DBG, FILE_name, "a");

		for (unsigned long jj = 0; jj < (sample_size_1ms); jj++)
		{
			fprintf(DataLog_DBG, "%lf ", Z_IFFT_DUMP[jj]);
		}
		fprintf(DataLog_DBG, "\n");
		fclose(DataLog_DBG);
#endif

		int a = 0;
		int* i = &a;
		// ��� �ִ밪 ���
		Peak1[f_CH.d_f_index] = getMaxdouble_test(Z_IFFT_DUMP, sample_size_1ms, i);

		Peak_index[f_CH.d_f_index] = a;

		// �ִ� ����� +- 10chip�� ���� 0�����ͷ� ��ȯ
		for (unsigned long jj = 0; jj < (sample_size_1ms); jj++)
		{
			if ((jj < (Peak_index[f_CH.d_f_index] + sample_size_1ms * 10.0 / 1023.0)) &&
				(jj > (Peak_index[f_CH.d_f_index] - sample_size_1ms * 10.0 / 1023.0)))
			{
				Z_IFFT_for_Peak2[jj] = 0;
			}
			else
			{
				Z_IFFT_for_Peak2[jj] = Z_IFFT_DUMP[jj];
			}
		}

		Peak2[f_CH.d_f_index] = getMaxdouble(Z_IFFT_for_Peak2, sample_size_1ms);

		Peak_Ratio[f_CH.d_f_index] = Peak1[f_CH.d_f_index] / Peak2[f_CH.d_f_index];

		f_CH.d_f_index++;

	}

	for (int i = 0; i < 10; i++)
	{
		MAX_INDEX = getMaxIndex(Peak_Ratio, 45);

		if (Peak_Ratio[MAX_INDEX] > FFT_acq_threshold)
		{
			f_CH.d_f_index = MAX_INDEX;
			f_CH.d_f = range_dopp_freq_GPS[f_CH.d_f_index];
			f_CH.s_c_FFT_peak_index = Peak_index[MAX_INDEX];
			f_CH.DCO_Carr_INC = f_CH.DCO_Carr_INC_IF + (int)((double)f_CH.d_f / DCO_Carr_RESOL);
			f_CH.s_c = round((1023.0 - (static_cast<double>(f_CH.s_c_FFT_peak_index) / static_cast<double>(sample_size_1ms) * 1023.0)) * 2) / 2;

			f_CH.codephase = f_CH.s_c_FFT_peak_index;
			f_CH.d_f_code = f_CH.d_f / 1540.0;

			// �Ӱ谪�� ������ ��� d_f, s_c�� ����ü�� ����
			f_CH.detect = true;
			f_CH.ACQ_Check = false;
			f_CH.TRC_Check = true;		////////////////////// ���� ///////////////////��ȣ ȹ�� off, ��ȣ ���� on
			detect_finish = 0;

			fopen_s(&DataLog_Acq, "DataLog_Acq_result.txt", "a");
			fprintf(DataLog_Acq, "SYS : GPS, PRN : %d, doppler freq : %f (Hz), code shift : %f (chip) \n", f_CH.PRN, f_CH.d_f, f_CH.s_c);
			fclose(DataLog_Acq);


			break;

		}
		else
		{
			Peak_Ratio[MAX_INDEX] = 0;
		}

	}

	if (detect_finish == 1)	//ȹ���� �������� ����
	{
		f_CH.d_f_index = 0;

		if (f_CH.PRN < 32) //MS // ���� ��ȣ ȹ�� ���� // �ش�ä�� �ʱ�ȭ �� - RPN ��� �籸��
		{
			f_CH.PRN = f_CH.PRN + 1;	//MS
			f_CH.DCO_Carr = 0;
			f_CH.DCO_Code = 0;
			f_CH.d_f = range_dopp_freq_GPS[0];
			f_CH.s_c = range_shift_code_GPS[0];
			f_CH.ACQ_Check = true;
			f_CH.Bit_Check = false;
		}
		else
		{
			f_CH.PRN = 0;
			f_CH.DCO_Carr = 0;
			f_CH.DCO_Code = 0;
			f_CH.d_f = range_dopp_freq_GPS[0];
			f_CH.s_c = range_shift_code_GPS[0];
			f_CH.ACQ_Check = false;
			f_CH.TRC_Check = false;
			f_CH.DCO_Code_Check = false;
			f_CH.Bit_Check = false;

		}
	}


	free(Z_IFFT_for_Peak2);
	free(Z_IFFT_DUMP);

#if ACQ_PLOT_2
	// Peak ratiio result
	fopen_s(&DataLog_Acq, "DataLog\\DataLog_Acq_Threshold.txt", "a");
	fprintf(DataLog_Acq, "%f\n", Peak_Ratio[MAX_INDEX]);
	fclose(DataLog_Acq);
#endif
}

void f_Peak_Detector_5ms(Channel_struct& f_CH, double* Z_IFFT)
{

	int target_PRN = 10;
	int target_d_f_index = 11;

	// 5ms peak detector(doppler search 100Hz)
	double Peak1[221] = { 0 };
	double Peak2[221] = { 0 };
	int Peak_index[221] = { 0 };
	double Peak_Ratio[221] = { 0 };
	int s_c_FFT_peak_index[221] = { 0 };
	double* Z_IFFT_for_Peak2 = NULL;
	Z_IFFT_for_Peak2 = (double*)malloc(sizeof(double) * (int)(f_s / 1000.0));
	double* Z_IFFT_DUMP;
	Z_IFFT_DUMP = (double*)malloc(sizeof(double) * (int)(f_s / 1000.0));
	int detect_finish = 1;
	int MAX_INDEX = 0;

	while (f_CH.d_f_index < 221)
	{

		for (unsigned long jj = 0; jj < (sample_size_1ms); jj++)
		{
			Z_IFFT_DUMP[jj] = Z_IFFT[((CorrSizeLocal * 5) * 221 * (f_CH.PRN - 1)) + (CorrSizeLocal * f_CH.d_f_index) + jj];
		}

#if ACQ_PLOT_1

		sprintf_s(FILE_name, "DataLog\\Datalog_Acq_Full_FFT_GPS%d.txt", f_CH.PRN);
		fopen_s(&DataLog_DBG, FILE_name, "a");

		for (unsigned long jj = 0; jj < (sample_size_1ms); jj++)
		{
			fprintf(DataLog_DBG, "%lf ", Z_IFFT_DUMP[jj]);
		}
		fprintf(DataLog_DBG, "\n");
		fclose(DataLog_DBG);
#endif

		int a = 0;
		int* i = &a;

		Peak1[f_CH.d_f_index] = getMaxdouble_test(Z_IFFT_DUMP, sample_size_1ms, i);

		Peak_index[f_CH.d_f_index] = a;

		// �ִ� ����� +- 10chip�� ���� 0�����ͷ� ��ȯ
		for (unsigned long jj = 0; jj < (sample_size_1ms); jj++)
		{
			if ((jj < (Peak_index[f_CH.d_f_index] + sample_size_1ms * 10.0 / 1023.0)) &&
				(jj > (Peak_index[f_CH.d_f_index] - sample_size_1ms * 10.0 / 1023.0)))
			{
				Z_IFFT_for_Peak2[jj] = 0;
			}
			else
			{
				Z_IFFT_for_Peak2[jj] = Z_IFFT_DUMP[jj];
			}
		}

		Peak2[f_CH.d_f_index] = getMaxdouble(Z_IFFT_for_Peak2, sample_size_1ms);

		Peak_Ratio[f_CH.d_f_index] = Peak1[f_CH.d_f_index] / Peak2[f_CH.d_f_index];

		f_CH.d_f_index++;

	}

	for (int i = 0; i < 10; i++)
	{
		MAX_INDEX = getMaxIndex(Peak_Ratio, 221);

		if (Peak_Ratio[MAX_INDEX] > FFT_acq_threshold)
		{
			f_CH.d_f_index = MAX_INDEX;
			f_CH.d_f = range_dopp_freq_GPS[f_CH.d_f_index];
			f_CH.s_c_FFT_peak_index = Peak_index[MAX_INDEX];
			f_CH.DCO_Carr_INC = f_CH.DCO_Carr_INC_IF + (int)((double)f_CH.d_f / DCO_Carr_RESOL);
			f_CH.s_c = round((1023.0 - (static_cast<double>(f_CH.s_c_FFT_peak_index) / static_cast<double>(sample_size_1ms) * 1023.0)) * 2) / 2;

			f_CH.codephase = f_CH.s_c_FFT_peak_index;
			f_CH.d_f_code = f_CH.d_f / 1540.0;

			// �Ӱ谪�� ������ ��� d_f, s_c�� ����ü�� ����
			f_CH.detect = true;
			f_CH.ACQ_Check = false;
			f_CH.TRC_Check = true;		////////////////////// ���� ///////////////////��ȣ ȹ�� off, ��ȣ ���� on
			detect_finish = 0;

			fopen_s(&DataLog_Acq, "DataLog_Acq_result.txt", "a");
			fprintf(DataLog_Acq, "SYS : GPS, PRN : %d, doppler freq : %f (Hz), code shift : %f (chip) \n", f_CH.PRN, f_CH.d_f, f_CH.s_c);
			fclose(DataLog_Acq);


			break;

		}
		else
		{
			Peak_Ratio[MAX_INDEX] = 0;
		}

	}

	if (detect_finish == 1)	//ȹ���� �������� ����
	{
		f_CH.d_f_index = 0;

		if (f_CH.PRN < 32) //MS // ���� ��ȣ ȹ�� ���� // �ش�ä�� �ʱ�ȭ �� - RPN ��� �籸��
		{
			f_CH.PRN = f_CH.PRN + 1;	//MS
			f_CH.DCO_Carr = 0;
			f_CH.DCO_Code = 0;
			f_CH.d_f = range_dopp_freq_GPS[0];
			f_CH.s_c = range_shift_code_GPS[0];
			f_CH.ACQ_Check = true;
			f_CH.Bit_Check = false;
		}
		else
		{
			f_CH.PRN = 0;
			f_CH.DCO_Carr = 0;
			f_CH.DCO_Code = 0;
			f_CH.d_f = range_dopp_freq_GPS[0];
			f_CH.s_c = range_shift_code_GPS[0];
			f_CH.ACQ_Check = false;
			f_CH.TRC_Check = false;
			f_CH.DCO_Code_Check = false;
			f_CH.Bit_Check = false;

		}
	}


	free(Z_IFFT_for_Peak2);
	free(Z_IFFT_DUMP);

#if ACQ_PLOT_2
	// Peak ratiio result
	fopen_s(&DataLog_Acq, "DataLog\\DataLog_Acq_Threshold.txt", "a");
	fprintf(DataLog_Acq, "%f\n", Peak_Ratio[MAX_INDEX]);
	fclose(DataLog_Acq);
#endif
}

void f_loop_filter(Channel_struct& f_CH, int DCO_Carr_bit, int DCO_Code_bit, double DCO_Carr_RESOL, double DCO_Code_RESOL, int f_s)	// FLL, DLL, PLL ����
{
	// FLL-PLL
	// Stage 1
	f_CH.PLL_B0[1] = (f_CH.PLL_coeff1 * (0.001) * f_CH.PLL_Discri) + (f_CH.FLL_coeff1 * (0.001) * f_CH.FLL_Discri) + f_CH.PLL_B0[0];
	f_CH.PLL_B0[2] = (f_CH.PLL_B0[1] + f_CH.PLL_B0[0]) * (0.5);
	f_CH.PLL_B0[0] = f_CH.PLL_B0[1];

	f_CH.PLL_1to2 = f_CH.PLL_B0[2] + (f_CH.PLL_Discri * f_CH.PLL_coeff2) + (f_CH.FLL_Discri * f_CH.FLL_coeff2);

	// Stage 2
	f_CH.PLL_B1[1] = f_CH.PLL_1to2 * (0.001) + f_CH.PLL_B1[0];
	f_CH.PLL_B1[2] = (f_CH.PLL_B1[1] + f_CH.PLL_B1[0]) * (0.5);
	f_CH.PLL_B1[0] = f_CH.PLL_B1[1];

	f_CH.PLL_2to3 = f_CH.PLL_B1[2] + (f_CH.PLL_Discri * f_CH.PLL_coeff3);

	f_CH.Filtered_PLL_Discri = f_CH.PLL_2to3 * 1 / (2 * M_PI) * pow((double)2, (double)DCO_Carr_bit) / f_s;



	f_CH.DCO_Carr_INC = f_CH.DCO_Carr_INC_IF + (int)((double)f_CH.d_f / DCO_Carr_RESOL) + f_CH.Filtered_PLL_Discri;
	f_CH.DCO_Code_INC = (int)(f_CH.DCO_Code_INC_IF + ((f_CH.DCO_Carr_INC - f_CH.DCO_Carr_INC_IF) / 1540) + (f_CH.DLL_Discri * f_CH.DLL_w0 / DCO_Code_RESOL * 20) + 0.5);
	f_CH.Doppler = ((f_CH.DCO_Carr_INC - f_CH.DCO_Carr_INC_IF) * DCO_Carr_RESOL);


}

void f_Data_Read(int ifdata, int TRK_first, int ACTIVECHANNEL, short* TRK_CH_IF_DATA, short* TRK_IN_IF_DATA, int* N_trk, bool* Trk_ready, int* remaining_data)
{
	if (TRK_first == 1) //ù TRK ���� �� 2�ֱ� IF DATA �о��
	{


		for (int i = 0; i < ACTIVECHANNEL; i++)
		{
			TRK_CH_IF_DATA[i * (2 * CorrSizeLocal) + (*N_trk)] = ifdata;

		}

		(*N_trk)++;

		if ((*N_trk) == (CorrSizeLocal * 2))
		{
			(*Trk_ready) = true;
			(*N_trk) = 0;
		}


	}
	else  //ù��° TRK �Ϸ� �� 1�ֱ� IF DATA �о��
	{
		TRK_IN_IF_DATA[(*remaining_data) + (*N_trk)] = ifdata;

		(*N_trk)++;

		if ((*N_trk) == (CorrSizeLocal))
		{
			(*Trk_ready) = true;
			(*remaining_data) = (*remaining_data) + (*N_trk);	//buffer�� �����ִ� ������ ����
			(*N_trk) = 0;
		}

		if ((*Trk_ready) == true)
		{

			int size_min = CorrSizeLocal;

			for (int i = 0; i < ACTIVECHANNEL; i++)
			{
				CH[i].valid_data = CH[i].valid_data - (CH[i].Second_2046 + 1);	//�� ä�� ���ۿ� ���� ������ ����

				memmove(TRK_CH_IF_DATA + (i * CorrSizeLocal * 2), TRK_CH_IF_DATA + (i * CorrSizeLocal * 2) + (CH[i].Second_2046 + 1), sizeof(short) * CH[i].valid_data);	//�� ä�� ���ۺ��� ������� ���� �����͸� ������ �̵�


				if (CH[i].valid_data < CorrSizeLocal)		//�� ���ۿ� ���� �����Ͱ� 1ms ���� ����(25000)��ŭ ���� ������ ������ ���� ���ۿ��� ������ ���
				{
					memcpy(TRK_CH_IF_DATA + (i * CorrSizeLocal * 2) + CH[i].valid_data, TRK_IN_IF_DATA + CH[i].used_data, sizeof(short) * (CorrSizeLocal - CH[i].valid_data));

					CH[i].used_data = CH[i].used_data + (CorrSizeLocal - CH[i].valid_data);
					CH[i].valid_data = CorrSizeLocal;

				}

				if (CH[i].used_data < size_min)
				{
					size_min = CH[i].used_data;
				}


			}

			(*remaining_data) = (*remaining_data) - size_min;


			memmove(TRK_IN_IF_DATA, TRK_IN_IF_DATA + size_min, sizeof(short) * (*remaining_data));		//������ ���ۿ� ����� �����͸� ������ ���� �����͸� ������ �̵�


			for (int i = 0; i < ACTIVECHANNEL; i++)
			{
				CH[i].used_data = CH[i].used_data - size_min;

			}

		}

	}

}

void f_1ms_Save(Channel_struct& f_CH)
{

	f_CH.bf_carr = f_CH.DCO_Carr;
	f_CH.bf_code = f_CH.DCO_Code;
	f_CH.bf_code_count = f_CH.Code_Count;
	f_CH.bf_ephoch_count = f_CH.Epoch_Count;
	f_CH.bf_ephoch_count_20 = f_CH.Epoch_Count_20m;
	f_CH.bf_code_index = f_CH.code_index;
}

void f_1ms_Read(Channel_struct& f_CH)
{

	f_CH.DCO_Carr = f_CH.bf_carr;
	f_CH.DCO_Code = f_CH.bf_code;
	f_CH.Code_Count = f_CH.bf_code_count;
	f_CH.Epoch_Count = f_CH.bf_ephoch_count;
	f_CH.Epoch_Count_20m = f_CH.bf_ephoch_count_20;
	f_CH.code_index = f_CH.bf_code_index;

}

void f_Code_Count(Channel_struct& f_CH, int TRK_first, int* first_check, int lnN, int* check_2046, int* remaining_data, short* Cur_IF_2bytes_sign, FILE* pFile_IFData, short* TRK_IN_IF_DATA, short* TRK_CH_IF_DATA, int k)
{
	if (TRK_first == 1)
	{
		if (f_CH.DCO_Code_Check == true)
		{
			f_CH.DCO_Code_Check_TRC = true;

			f_CH.code_index += 1;
			f_CH.DCO_Code_Check = false;


			f_CH.Code_Count++;	// 1Chip ���� ī����

			if (f_CH.Code_Count == 2046) // 1023chip �������� ������Ʈ
			{

				f_CH.Epoch_Count++;
				f_CH.Code_Count = 0;
				if (f_CH.Epoch_Count == 20) //*20ms �������� ������Ʈ
				{
					f_CH.Epoch_Count_20m++;
					f_CH.Epoch_Count = 0;

					if (f_CH.Epoch_Count_20m == 50) //*1s�϶� �ʱ�ȭ
						f_CH.Epoch_Count_20m = 0;
				}

			}

			if ((int)(f_CH.code_index + f_CH.s_c * 2) % 2046 == 0 && (*first_check) == 0)
			{
				f_CH.First_2046 = lnN;

				(*first_check)++;
			}
			else if (f_CH.code_index == (2046) && (*first_check) == 2)
			{
				f_CH.Second_2046 = lnN;

				f_CH.corr_size = CorrSizeLocal;		//corr ���� 2ms ������ �� 1ms�� ���� �����̱� ������ CorrSizeLocal�� ����
				f_CH.valid_data = CorrSizeLocal;
				(*first_check)++;

				f_1ms_Save(f_CH);

				f_CH.onems_count += lnN + 1;		//�ι�° 1ms �����Ϳ��� ����� �����Ͱ� ���������� Ȯ��
			}

		}
		else
		{
			f_CH.DCO_Code_Check_TRC = false;
		}
	}
	else
	{
		if (f_CH.DCO_Code_Check == true)
		{
			f_CH.DCO_Code_Check_TRC = true;

			f_CH.code_index += 1;
			f_CH.DCO_Code_Check = false;

			f_CH.Code_Count++;	// 1Chip ���� ī����

			if (f_CH.Code_Count == 2046) // 1023chip �������� ������Ʈ
			{
				f_CH.Epoch_Count++;
				f_CH.Code_Count = 0;
				if (f_CH.Epoch_Count == 20) //*20ms �������� ������Ʈ
				{
					f_CH.Epoch_Count_20m++;
					f_CH.Epoch_Count = 0;

					if (f_CH.Epoch_Count_20m == 50) //*1s�϶� �ʱ�ȭ
						f_CH.Epoch_Count_20m = 0;
				}


			}

			if (f_CH.code_index == (2046))
			{
				int ifdata = 0;

				f_CH.First_2046 = 0;
				f_CH.Second_2046 = lnN;

				f_CH.corr_size = lnN + 1;
				(*check_2046) = 1;
				f_CH.code_index_ready = 1;

				f_1ms_Save(f_CH);

				if (f_CH.valid_data < f_CH.corr_size)		//�� ä�� ���ۿ� �ִ� �����Ͱ� corr_size���� ������ ������ ����
				{

					if (((*remaining_data) - f_CH.used_data) < (f_CH.corr_size - f_CH.valid_data))
					{
						for (int i = 0; i < CorrSizeLocal; i++)	//���ۿ� 1ms ������ �о��
						{
							fread(Cur_IF_2bytes_sign, 1, sizeof((*Cur_IF_2bytes_sign)), pFile_IFData);
							ifdata = (*Cur_IF_2bytes_sign);
							TRK_IN_IF_DATA[(*remaining_data) + i] = ifdata;
						}

						(*remaining_data) = (*remaining_data) + CorrSizeLocal;

					}
					memcpy(TRK_CH_IF_DATA + (k * CorrSizeLocal * 2) + CorrSizeLocal, TRK_IN_IF_DATA + f_CH.used_data, sizeof(short) * (f_CH.corr_size - f_CH.valid_data));

					f_CH.used_data = f_CH.used_data + (f_CH.corr_size - f_CH.valid_data);
					f_CH.valid_data = f_CH.valid_data + (f_CH.corr_size - f_CH.valid_data);
				}


			}


		}
		else
		{
			CH[k].DCO_Code_Check_TRC = false;// DCO_Code_Check on ������ �ƴѰ�� false�� ����																			
		}
	}
}

void f_BitCheck(Channel_struct& f_CH, double temp_phase_IQ)
{
	int Bit_Change_threshold = 20;
	int Bit_Change_search_index = 0;

	//��Ʈ���Ⱑ �Ϸ��� ����
	if (f_CH.Bit_Check == true)
	{
		if (f_CH.Bit_Change_Count_index_20msec == f_CH.Bit_Change_Location)
		{
			if (temp_phase_IQ > M_PI_2 && temp_phase_IQ < (3 * M_PI_2))
			{
				f_CH.Bit_Change_Count[f_CH.Bit_Change_Location]++;
				f_CH.Navi_Data_now = -1 * (f_CH.Navi_Data_pre - 1);
			}
			else {
				f_CH.Navi_Data_now = f_CH.Navi_Data_pre;
			}

			f_CH.Navi_Data_pre = f_CH.Navi_Data_now;

			f_CH.Navi_Data_30bit[0] = ((f_CH.Navi_Data_30bit[0] << 1) | (Data_Emit_32bit(f_CH.Navi_Data_30bit[1], 1, 1) & 0x01));
			f_CH.Navi_Data_30bit[1] = ((f_CH.Navi_Data_30bit[1] << 1) | f_CH.Navi_Data_now);

			// ������ ���� �Ϸ�
			if (f_CH.Frame_Sync_Check == true)
			{
				if (f_CH.Navi_Data_index == 29)
				{

					f_CH.Navi_Data_index = -1;
					f_CH.GPS_word[f_CH.GPS_word_index] = f_CH.Navi_Data_30bit[1];

					// ���κ�Ʈ ���� Ȯ��
					if (f_CH.GPS_word[f_CH.GPS_word_index] & 0x40000000L)
					{
						f_CH.GPS_word[f_CH.GPS_word_index] ^= 0x3FFFFFC0L;
					}

					// �и�Ƽ üũ 
					f_Parity_Check(f_CH.GPS_word[f_CH.GPS_word_index], f_CH.Parity_Check);

					if (f_CH.Parity_Check == false)
					{
						f_CH.Parity_Fail_Count++;
					}

					// Parity Check && Monitoring
					if (f_CH.Parity_Fail_Count >= 1)
					{
						f_CH.Parity_Fail_Count = 0;
						f_CH.Frame_Sync_Check = false;
						f_CH.Preamble_Check = false;
						f_CH.Bit_Inversion = false;
						f_CH.Zero_bits_Check = false;
						f_CH.Parity_Check = false;
						f_CH.Parity_Fail_Count = 0;
						f_CH.Frame_init_pos_Check = false;
					}
					// ���� word index ������Ʈ
					if (f_CH.GPS_word_index == 9)
					{
						f_CH.Subframe_ID_Check = Data_Emit_32bit(f_CH.GPS_word[1], 20, 3);

						// ������ ���ڵ� ����
						f_Data_Decoding(f_CH);
						f_CH.GPS_word_index = 0;
					}
					else {
						f_CH.GPS_word_index++;
					}
				}
			}
			// ������ ���� �̿Ϸ�
			else // <------������ ��ũ X
			{
				if (f_CH.Navi_Data_index == 59)	 // 60 Bit ������ ���� �� ���
				{
					// ��ȿ�� �˻� ���� 
					f_GPS_Word0_Word1_check(f_CH.Navi_Data_30bit[0], f_CH.Navi_Data_30bit[1], f_CH.Frame_Sync_Check, f_CH);

					// ������ ���Ⱑ �Ϸ�Ǿ��ٸ� 
					if (f_CH.Frame_Sync_Check == true)
					{
						// Word0, Word1�� Pairty Decoding (D1, D24 bit�� ���� word�� D30 Ȯ�� �� ����)
						if (f_CH.Navi_Data_30bit[0] & 0x40000000L)
							f_CH.Navi_Data_30bit[0] ^= 0x3FFFFFC0L;
						if (f_CH.Navi_Data_30bit[1] & 0x40000000L)
							f_CH.Navi_Data_30bit[1] ^= 0x3FFFFFC0L;

						// GPS word�� ������ ����
						f_CH.GPS_word[0] = f_CH.Navi_Data_30bit[0];
						f_CH.GPS_word[1] = f_CH.Navi_Data_30bit[1];

						// GPS word index ���� 
						f_CH.GPS_word_index = 2;
						f_CH.Navi_Data_index = -1;

						if (f_CH.Frame_init_pos_Check == false)
						{
							f_CH.Frame_init_pos = 49 - ((f_CH.Epoch_Count_20m + 1 + 40) % 50);

							if (f_CH.Bit_init_pos == 0)
							{
								f_CH.Frame_init_pos++;
								if (f_CH.Frame_init_pos > 49)
									f_CH.Frame_init_pos -= 50;
							}

							f_CH.Frame_init_pos_Check = true;
						}



					}
					else
					{
						f_CH.Frame_Sync_Check = false;
						f_CH.Preamble_Check = false;
						f_CH.Bit_Inversion = false;
						f_CH.Zero_bits_Check = false;
						f_CH.Parity_Check = false;
						f_CH.Frame_init_pos_Check = false;
						f_CH.Navi_Data_index--;
					}
				}
			}
			f_CH.Navi_Data_index++;
		}


		// �޽��� 1��Ʈ(20ms) ������ ī����
		if (f_CH.Bit_Change_Count_index_20msec == 19)
		{
			f_CH.Bit_Change_Count_index_20msec = 0;

			// 5�ʸ��� Bit_Chagne_Count ����͸� �˻�
			if (f_CH.Bit_Change_Count_index_1000msec == 50 * 5 - 1)
			{
				f_CH.Bit_Change_Count_index_1000msec = 0;


				if (f_CH.Bit_Change_Count[f_CH.Bit_Change_Location] < Bit_Change_threshold)
				{
					f_CH.Bit_Check = false;
					f_CH.Frame_init_pos_Check = false;
				}

				f_CH.Bit_Change_Count[f_CH.Bit_Change_Location] = 0;
			}
			else
			{
				f_CH.Bit_Change_Count_index_1000msec++;
			}
		}
		else
		{
			f_CH.Bit_Change_Count_index_20msec++;
		}
	}


	// ��Ʈ ��ũ �̰��� ����
	//if (CH[k].Bit_Check == false)
	else
	{
		if (temp_phase_IQ > M_PI_2 && temp_phase_IQ < (3 * M_PI_2))
		{
			//	�޽��� 1��Ʈ ������ index(0~19)�� Bit ��ȯ ī��Ʈ
			f_CH.Bit_Change_Count[f_CH.Bit_Change_Count_index_20msec]++;
		}

		if (f_CH.Bit_Change_Count_index_20msec == 19)
		{
			f_CH.Bit_Change_Count_index_20msec = 0;

			if (f_CH.Bit_Change_Count_index_1000msec == 49)
			{
				f_CH.Bit_Change_Count_index_1000msec = 0;

				for (Bit_Change_search_index = 0; Bit_Change_search_index < 20; Bit_Change_search_index++)
				{
					if (f_CH.Bit_Change_Count[Bit_Change_search_index] > Bit_Change_threshold)
					{
						//fopen_s(&DataLog_DBG, "DataLog\\DataLog_bit_change_check.txt", "a");
						//fprintf(DataLog_DBG, "GPS[%d] -> Bit_Change_index[%d] = %d\n",
						//	f_CH.PRN, Bit_Change_search_index, f_CH.Bit_Change_Count[Bit_Change_search_index]);
						//fclose(DataLog_DBG);

						f_CH.Bit_Change_Location = Bit_Change_search_index;
						f_CH.Bit_Check = true;
						f_CH.Bit_init_pos = 19 - f_CH.Bit_Change_Location;

						break;
					}
				}
			}
			else
			{
				f_CH.Bit_Change_Count_index_1000msec++;
			}
		}
		else
		{
			f_CH.Bit_Change_Count_index_20msec++;
		}
	}
}

void f_Tracking_GPS(int ifdata, int* TRK_first, int ACTIVECHANNEL, short* TRK_CH_IF_DATA, short* TRK_IN_IF_DATA, int* N_trk, bool* Trk_ready, int* remaining_data, short* Cur_IF_2bytes_sign, FILE* pFile_IFData, int* Gpu_CA_code_GPS_2046, int* CorrSizeLocal_count)
{

	int k = 0;
	int DCO_Carr_temp = 0;
	//int CorrSizeLocal_count = 0;
	int trk_blocks = 0;
	double temp_phase_IQ = 0.0;
	double dot = 0;
	double cross = 0;

	f_Data_Read(ifdata, (*TRK_first), ACTIVECHANNEL, TRK_CH_IF_DATA, TRK_IN_IF_DATA, N_trk, Trk_ready, remaining_data);


	if ((*Trk_ready) == true)
	{

		for (k = 0; k < ACTIVECHANNEL; k++) //���� ä�� 7�� ����
		{

			if (CH[k].TRC_Check == true)
			{
				Cpu_trk_PRN[k] = (short)CH[k].PRN;

				int check_2046 = 0;

				//int first_check = 0;

				if ((*TRK_first) == 1)
				{
					int first_check = 0;

					for (int lnN = 0; lnN < CorrSizeLocal * 2; lnN++)
					{

						f_Carr_DCO(CH[k].DCO_Carr, g_DCO_Carr_bit, CH[k].DCO_Carr_INC);
						DCO_Carr_temp = (0x00000007 & (CH[k].DCO_Carr >> (g_DCO_Carr_bit - 3)));
						f_Code_DCO(CH[k].DCO_Code, g_DCO_Code_bit, CH[k].DCO_Code_INC * 2, CH[k].DCO_Code_Check, CH[k].PRN);

						f_Code_Count(CH[k], (*TRK_first), &first_check, lnN, &check_2046, remaining_data, Cur_IF_2bytes_sign, pFile_IFData, TRK_IN_IF_DATA, TRK_CH_IF_DATA, k);

						Cpu_carr_temp[(k * CorrSizeLocal * 2) + lnN] = (short)DCO_Carr_temp;
						Cpu_code_index[(k * CorrSizeLocal * 2) + lnN] = (short)CH[k].code_index;

						if (first_check == 1 || first_check == 3)
						{
							CH[k].code_index = 0;
							first_check++;
						}

					}	// <------- for (int lnN = 0; lnN < CorrSizeLocal * 2; lnN++)

					for (int i = 0; i < ACTIVECHANNEL; i++)
					{
						f_1ms_Read(CH[i]);
					}
				}
				else     //if (TRK_first == 1)�� else
				{
					int replica_gen_check = 0;
					int lnN = 0;
					check_2046 = 0;

					f_1ms_Read(CH[k]);

					CH[k].code_index = 0;


					for (lnN = 0; lnN < CorrSizeLocal + 100; lnN++)		//���÷� ȿ���� corr size�� �þ �� �ֱ� ������ 1ms ���̿� 100������ ������ ��
					{
						if (check_2046 == 1)	//code ���� 1ms�� �Ǹ� ���Ĵ� 0���� ����
						{
							CH[k].code_index = 0;
							DCO_Carr_temp = 0;
						}
						else
						{
							f_Carr_DCO(CH[k].DCO_Carr, g_DCO_Carr_bit, CH[k].DCO_Carr_INC);
							DCO_Carr_temp = (0x00000007 & (CH[k].DCO_Carr >> (g_DCO_Carr_bit - 3)));

							f_Code_DCO(CH[k].DCO_Code, g_DCO_Code_bit, CH[k].DCO_Code_INC * 2, CH[k].DCO_Code_Check, CH[k].PRN);

							f_Code_Count(CH[k], (*TRK_first), 0, lnN, &check_2046, remaining_data, Cur_IF_2bytes_sign, pFile_IFData, TRK_IN_IF_DATA, TRK_CH_IF_DATA, k);

							CH[k].onems_count++;

						}

						Cpu_carr_temp[(k * CorrSizeLocal * 2) + lnN] = (short)DCO_Carr_temp;
						Cpu_code_index[(k * CorrSizeLocal * 2) + lnN] = (short)CH[k].code_index;

						if (CH[k].onems_count == NAV.rec_TIC_Count_Interval + 1 && CH[k].code_index != (2046))
						{
							CH[k].onems_ready = 1;

							memcpy(&NAV_CH[k], &CH[k], sizeof(struct Channel_struct));

							CH[k].onems_count = 0;

							(*CorrSizeLocal_count) += 1;
							//printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
							//printf("trk_test = %d rec_TIC_Count = %d\n", trk_test, CorrSizeLocal_count);

						}
					}

				}

			}
		}

		for (int i = 0; i < ACTIVECHANNEL; i++)
		{
			Cpu_trk_CH_2046[2 * i] = CH[i].First_2046;
			Cpu_trk_CH_2046[2 * i + 1] = CH[i].Second_2046;
			Cpu_s_c[i] = CH[i].s_c * 2;

			f_1ms_Read(CH[i]);

		}

		cudaMemcpy(Gpu_trk_ch_ifdata, TRK_CH_IF_DATA, sizeof(short) * CorrSizeLocal * 2 * ACTIVECHANNEL, cudaMemcpyHostToDevice);
		cudaMemcpy(Gpu_trk_PRN, Cpu_trk_PRN, sizeof(short) * ACTIVECHANNEL, cudaMemcpyHostToDevice);
		cudaMemcpy(Gpu_code_index, Cpu_code_index, sizeof(short) * (CorrSizeLocal * 2 * ACTIVECHANNEL), cudaMemcpyHostToDevice);
		cudaMemcpy(Gpu_carr_temp, Cpu_carr_temp, sizeof(short) * (CorrSizeLocal * 2 * ACTIVECHANNEL), cudaMemcpyHostToDevice);
		cudaMemcpy(Gpu_s_c, Cpu_s_c, sizeof(int) * ACTIVECHANNEL, cudaMemcpyHostToDevice);
		cudaMemcpy(Gpu_trk_CH_2046, Cpu_trk_CH_2046, sizeof(int) * ACTIVECHANNEL * 2, cudaMemcpyHostToDevice);

		//GPU ������ �޸�
		//cudaMemcpy(Gpu_trk_count, Cpu_trk_count, sizeof(int) * 1, cudaMemcpyHostToDevice);

		//cuda_op.f_Correlation // GPU ��� ���� �Լ� ȣ��
		f_Correlation // GPU ��� ���� �Լ� ȣ��
		(
			Gpu_trk_PRN,
			Gpu_code_index,
			Gpu_carr_temp,
			Gpu_s_c,
			Gpu_trk_ch_ifdata,
			Gpu_out_partial_dump,
			Gpu_CA_code_GPS_2046,
			ACTIVECHANNEL,
			Gpu_trk_CH_2046,
			(*TRK_first),
			Gpu_trk_count

		);
		(*Cpu_trk_count)++;

		if ((*TRK_first) == 1)
		{
			cudaMemcpy(Cpu_out_partial_dump, Gpu_out_partial_dump, sizeof(float) * First_trk_blocksPerGrid * 6 * ACTIVECHANNEL, cudaMemcpyDeviceToHost);
			trk_blocks = First_trk_blocksPerGrid;


			for (int i = 0; i < ACTIVECHANNEL; i++)
			{
				memmove(TRK_CH_IF_DATA + (i * CorrSizeLocal * 2), TRK_CH_IF_DATA + (i * CorrSizeLocal * 2) + CorrSizeLocal, sizeof(short) * CorrSizeLocal);
				CH[i].Second_2046 = CH[i].Second_2046 - CorrSizeLocal;

			}

		}
		else
		{
			cudaMemcpy(Cpu_out_partial_dump, Gpu_out_partial_dump, sizeof(float) * trk_blocksPerGrid * 6 * ACTIVECHANNEL, cudaMemcpyDeviceToHost);
			trk_blocks = trk_blocksPerGrid;
		}

		

		for (k = 0; k < ACTIVECHANNEL; k++) //���� ä�� 7�� ����
		{

			if (CH[k].TRC_Check == true)
			{
	
				CH[k].temp_dump_I_E = 0;
				CH[k].temp_dump_I_P = 0;
				CH[k].temp_dump_I_L = 0;
				CH[k].temp_dump_Q_E = 0;
				CH[k].temp_dump_Q_P = 0;
				CH[k].temp_dump_Q_L = 0;

				for (int iblock = 0; iblock < trk_blocks; iblock++)
				{
					CH[k].temp_dump_I_E += Cpu_out_partial_dump[(k * trk_blocks * 6) + iblock];
					CH[k].temp_dump_I_P += Cpu_out_partial_dump[(k * trk_blocks * 6) + (trk_blocks * 1) + iblock];
					CH[k].temp_dump_I_L += Cpu_out_partial_dump[(k * trk_blocks * 6) + (trk_blocks * 2) + iblock];
					CH[k].temp_dump_Q_E += Cpu_out_partial_dump[(k * trk_blocks * 6) + (trk_blocks * 3) + iblock];
					CH[k].temp_dump_Q_P += Cpu_out_partial_dump[(k * trk_blocks * 6) + (trk_blocks * 4) + iblock];
					CH[k].temp_dump_Q_L += Cpu_out_partial_dump[(k * trk_blocks * 6) + (trk_blocks * 5) + iblock];

				}
				CH[k].dump_finish = 1;
				

				if (CH[k].dump_finish == 1)
				{
					if (CH[k].Bit_Check == true)
					{
						CN0_Estimate(CH[k], CH[k].temp_dump_I_P, CH[k].temp_dump_Q_P, M_max, K_max);
					}

					//(*TRK_first) = 0;

					CH[k].trk_2046_count++;

					

					CH[k].DCO_Code_Check_TRC = false;
					CH[k].code_index = CH[k].bf_code_index - 2046;
					CH[k].dump_finish = 0;
					(*TRK_first) = 0;

					CH[k].DCO_Carr = CH[k].bf_carr;
					CH[k].DCO_Code = CH[k].bf_code;
					CH[k].Code_Count = 0;
					CH[k].Epoch_Count = CH[k].bf_ephoch_count;
					CH[k].Epoch_Count_20m = CH[k].bf_ephoch_count_20;


					CH[k].Z = (int)((CH[k].temp_dump_I_P * CH[k].temp_dump_I_P + CH[k].temp_dump_Q_P * CH[k].temp_dump_Q_P) / (f_s / 1000));

					CH[k].DLL_Discri = (((double)sqrt(CH[k].temp_dump_I_E * CH[k].temp_dump_I_E + CH[k].temp_dump_Q_E * CH[k].temp_dump_Q_E)
						- (double)sqrt(CH[k].temp_dump_I_L * CH[k].temp_dump_I_L + CH[k].temp_dump_Q_L * CH[k].temp_dump_Q_L))
						/ ((double)sqrt(CH[k].temp_dump_I_E * CH[k].temp_dump_I_E + CH[k].temp_dump_Q_E * CH[k].temp_dump_Q_E)
							+ (double)sqrt(CH[k].temp_dump_I_L * CH[k].temp_dump_I_L + CH[k].temp_dump_Q_L * CH[k].temp_dump_Q_L))) / 2.0;

					CH[k].PLL_Discri = atan((double)CH[k].temp_dump_Q_P / (double)CH[k].temp_dump_I_P);

					temp_phase_IQ = abs(atan2((double)CH[k].temp_dump_Q_P, (double)CH[k].temp_dump_I_P) - atan2((double)CH[k].temp_dump_Q_P_pre, (double)CH[k].temp_dump_I_P_pre));


					if (CH[k].temp_dump_I_P_pre != 0 && CH[k].temp_dump_Q_P_pre != 0)
					{
						if ((temp_phase_IQ > M_PI_2 && temp_phase_IQ < (3 * M_PI_2)) && CH[k].FLL_Discri_inv_start == true)
						{
							dot = (double)(-1 * CH[k].temp_dump_I_P * CH[k].temp_dump_I_P_pre - CH[k].temp_dump_Q_P * CH[k].temp_dump_Q_P_pre);
							cross = (double)(-1 * CH[k].temp_dump_I_P_pre * CH[k].temp_dump_Q_P + CH[k].temp_dump_Q_P_pre * CH[k].temp_dump_I_P);
						}

						else
						{
							dot = (double)(CH[k].temp_dump_I_P * CH[k].temp_dump_I_P_pre + CH[k].temp_dump_Q_P * CH[k].temp_dump_Q_P_pre);
							cross = (double)(CH[k].temp_dump_I_P_pre * CH[k].temp_dump_Q_P - CH[k].temp_dump_Q_P_pre * CH[k].temp_dump_I_P);
						}

						if (CH[k].FLL_Discri_Count > 100)
						{
							CH[k].FLL_Discri_inv_start = true;
						}
						else {
							CH[k].FLL_Discri_Count++;
						}

						CH[k].FLL_Discri = atan2(cross, dot) / (0.001) / (2 * 3.14);
					}

					CH[k].Z_5 += CH[k].Z;
					CH[k].Z_5_index++;
					if (CH[k].Z_5_index == 5)
					{

						if (CH[k].Z_5 / 5 < (Z_threshold / 5))
						{

							//fopen_s(&DataLog_Acq, "DataLog\\DataLog_Acq.txt", "a");
							//fprintf(DataLog_Acq, "TRK_Fail SYS : GPS, PRN : %d, doppler freq : %f (Hz), code shift : %f (chip) \n", CH[k].PRN, CH[k].d_f, CH[k].s_c);
							//fclose(DataLog_Acq);

						}
						CH[k].Z_5 = 0;
						CH[k].Z_5_index = 0;
					}

					f_loop_filter(CH[k], g_DCO_Carr_bit, g_DCO_Code_bit, DCO_Carr_RESOL, DCO_Code_RESOL, f_s);

					f_BitCheck(CH[k], temp_phase_IQ);

#if TRK_PLOT	

					sprintf_s(FILE_name, "DataLog\\DataLog_Trk_GPS%d.txt", CH[k].PRN);
					fopen_s(&DataLog_DBG, FILE_name, "a");
					fprintf(DataLog_DBG, "%lld %f %f %f %f %f %f %f %f %lld %lld %d\n",
						(long long)((CH[k].temp_dump_I_P * CH[k].temp_dump_I_P + CH[k].temp_dump_Q_P * CH[k].temp_dump_Q_P) / (f_s / 100000)),
						CH[k].DLL_Discri, CH[k].FLL_Discri,
						CH[k].Filtered_DLL_Discri, (CH[k].Filtered_FLL_Discri + CH[k].Filtered_PLL_Discri) * DCO_Carr_RESOL,
						CH[k].PLL_Discri, CH[k].PLL_A2[1], CH[k].temp_dump_I_P, CH[k].temp_dump_Q_P,
						(long long)((CH[k].temp_dump_I_E * CH[k].temp_dump_I_E + CH[k].temp_dump_Q_E * CH[k].temp_dump_Q_E) / (f_s / 100000)),
						(long long)((CH[k].temp_dump_I_L * CH[k].temp_dump_I_L + CH[k].temp_dump_Q_L * CH[k].temp_dump_Q_L) / (f_s / 100000)),
						CH[k].DCO_Code_INC
					);
					fclose(DataLog_DBG);
#endif

					CH[k].temp_dump_I_P_pre = CH[k].temp_dump_I_P;
					CH[k].temp_dump_Q_P_pre = CH[k].temp_dump_Q_P;

					CH[k].temp_dump_I_E = 0;
					CH[k].temp_dump_Q_E = 0;
					CH[k].temp_dump_I_P = 0;
					CH[k].temp_dump_Q_P = 0;
					CH[k].temp_dump_I_L = 0;
					CH[k].temp_dump_Q_L = 0;

				}

				if (CH[k].onems_count == NAV.rec_TIC_Count_Interval + 1 && CH[k].code_index_ready == 1)
				{
					memcpy(&NAV_CH[k], &CH[k], sizeof(struct Channel_struct));		
					(*CorrSizeLocal_count) += 1;
					CH[k].onems_count = 0;
					CH[k].code_index_ready = 0;
				}


			} 
		} 

		if ((*CorrSizeLocal_count) == ACTIVECHANNEL)
		{
			NAV.rec_TIC_Count += NAV.rec_TIC_Count_Interval;
			(*CorrSizeLocal_count) = 0;
		}



	} 

	(*Trk_ready) = false;

}