#pragma once

#include "Structs.h"
#include "Global.h"
#include "DSP.h"
#define _USE_MATH_DEFINES

#include "math.h"

void f_GNSS_SDR_Init(Channel_struct *f_CH, Navigation_struct &NAV);	 // 수신기 초기화
