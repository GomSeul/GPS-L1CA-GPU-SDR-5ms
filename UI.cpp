#include "UI.h"
#include "stdafx.h"
#include "Structs.h"
#include "Global.h"
#define _USE_MATH_DEFINES
#include "math.h"

void f_Console_Plot_GPS(Channel_struct* CH, Navigation_struct& NAV)
{

	if(ui_check == 1)
	{
		ui_check = 0;
		double update_time = (double)NAV.rec_TIC_Count / (double)f_s;

		printf("\n\n");
		printf("\x1B[1m[TIC:%9.3lf] [GPS_time:%9.3lf] [Process : %3.2f%%]\x1B[0m\n", NAV.rec_TIC + update_time, NAV.GPS_rec_time + update_time, ((NAV.rec_TIC + update_time) / SIGNAL_Time) * 100);
		printf("===========================================================================\n");
		printf("\x1B[1mECEF       :\x1B[0m %15.6f[m]   %15.6f[m]   %15.6f[m]\n", NAV.pos_ecef.x, NAV.pos_ecef.y, NAV.pos_ecef.z);
		printf("\x1B[1mLLH        :\x1B[0m %15.6f[deg] %15.6f[deg] %15.6f[m]\n", NAV.pos_llh.lat * (180.0 / M_PI), NAV.pos_llh.lon * (180.0 / M_PI), NAV.pos_llh.hgt);
		printf("\x1B[1mM_LLH      :\x1B[0m %15.6f[deg] %15.6f[deg] %15.6f[m]\n", (NAV.Mean_pos_llh.lat * (180.0 / M_PI)) / NAV_count, (NAV.Mean_pos_llh.lon * (180.0 / M_PI)) / NAV_count, NAV.Mean_pos_llh.hgt / NAV_count);
		printf("\x1B[1mVel[NED]   :\x1B[0m %15.6f[m/s] %15.6f[m/s] %15.6f[m/s]\n", NAV.vel_neu.x, NAV.vel_neu.y, NAV.vel_neu.z);
		printf("\x1B[1mVisible SV :\x1B[0m GPS %2d [#], SUM %2d [#]\n", NAV.GPS_L1_CA_NVS, NAV.GPS_L1_CA_NVS);
		printf("===========================================================================\n");
		printf("CH[#]   SYS   PRN   Elv[deg]   Azi[deg]   Dopp[Hz]   NCO[Hz]   SFID     PRerr[m]    LOCK   C/N0    |       X              Y              Z      |    Vel_X      Vel_Y      Vel_Z\n");
		printf("==================================================================================================================================================================================\n");

		for (int k = 0; k < Ch_num; k++) {
			if (CH[k].PRN > 0 && (CH[k].SYS == 0 || CH[k].SYS == 2))
			{
				printf("%3d", k);
				printf("%7d", CH[k].SYS);
				printf("%7d", CH[k].PRN);
				printf("%9.2f", CH[k].EL);
				printf("%11.2f", CH[k].AZ);
				printf("%12.2f", CH[k].d_f);
				printf("%11.2f", CH[k].Doppler);
				if (CH[k].SYS == 0)
					printf("%6d", CH[k].Subframe_ID_Check);
				else
					printf("%6d", 0);

				printf("%14.6f", (CH[k].Pred_Range - CH[k].Code_Meas));

				printf("   ");
				if (CH[k].ACQ_Check == true)
					printf("A");
				else
					printf(" ");

				if (CH[k].TRC_Check == true)
					printf("T");
				else
					printf(" ");

				if (CH[k].SYS == 0)
				{
					if (CH[k].Bit_Check == true)
						printf("B");
					else
						printf(" ");
					if (CH[k].Frame_Sync_Check == true)
						printf("F");
					else
						printf(" ");
					if (CH[k].Subframe_123_Check == true)
						printf("E");
					else
						printf(" ");
				}
				else
				{
					printf("  ");
				}

				printf("%8.2f", CH[k].CN0);
				printf("   |");
				printf("%14.2f", CH[k].SV_Pos.x);
				printf("%14.2f", CH[k].SV_Pos.y);
				printf("%14.2f", CH[k].SV_Pos.z);
				printf("  |");
				printf("%10.2f", CH[k].SV_Vel.x);
				printf("%11.2f", CH[k].SV_Vel.y);
				printf("%11.2f\n", CH[k].SV_Vel.z);

			}

		}


		COORD pos = { 0,0 };
		SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pos);
	}
}

