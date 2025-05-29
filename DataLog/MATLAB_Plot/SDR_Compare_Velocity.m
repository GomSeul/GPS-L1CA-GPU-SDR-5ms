clc;
clear;
close all;

D2R = pi/180;
R2D = 180/pi;


%% REFERENCE 궤적 저장
TRA_bin = fopen('[TRA] STATIC.bin');
A = fread(TRA_bin,'double');
len = length(A)/16; 
ref_temp = zeros(len,16); 

% 0.1 단위로 레퍼런스 저장
for i = 1:len
    idx = (i-1)*16 + 1;
    if mod(A(idx),0.1) == 0  && A(idx) 
        ref_temp(i,:) = A(idx:idx+15);
    end
end

% 0인 행 제거
ref_temp(~any(ref_temp,2),:) = [];
ref_temp = round(ref_temp*10^7)/10^7; % round to 7 decimal places
fclose(TRA_bin);

REF_TRA = [ref_temp(:,1:7), ref_temp(:,11:13)];

%% REFERENCE 궤적 출력

plot_3d_reftrj(REF_TRA(:,2:4),1,'m','m');
plot_reftrj(1,REF_TRA);

%% VELOCITY 오차 RMS 계산

% 1 : time
% 2~4 : LLH
% 5~7 : Velocity_NED

REF_VELOCITY = [ref_temp(:,1), ref_temp(:,5:7)];

fileID = fopen('../Datalog_Vel_Out.txt','r');
SDR_VEL_TEMP = textscan(fileID,'%f %f %f %f','Delimiter','\t');
fclose(fileID);
SDR_VEL_TEMP = cell2mat(SDR_VEL_TEMP);
SDR_VEL_TEMP = round(SDR_VEL_TEMP*10^7)/10^7; % round to 7 decimal places

MATCH = [];

for row_ref = 1:size(REF_VELOCITY, 1)
    ref_value_rounded = round(REF_VELOCITY(row_ref, 1), 1); % 첫 번째 자리까지 반올림
    matching_indices = find(ismember(round(SDR_VEL_TEMP(:,1), 1), ref_value_rounded));
    if ~isempty(matching_indices)
        MATCH = [MATCH, row_ref];
    end
end

REF_VELOCITY = REF_VELOCITY(MATCH, :);

VEL_ERR = REF_VELOCITY - SDR_VEL_TEMP;

TIME = SDR_VEL_TEMP(:, 1); 
VEL_ERR_N = VEL_ERR(:, 2); 
VEL_ERR_E = VEL_ERR(:, 3); 
VEL_ERR_D = VEL_ERR(:, 4); 

MEAN_VEL_ERR_N = mean(VEL_ERR_N); 
MEAN_VEL_ERR_E = mean(VEL_ERR_E); 
MEAN_VEL_ERR_D = mean(VEL_ERR_D); 

RMS_VEL_ERR_N = rms(VEL_ERR_N);
RMS_VEL_ERR_E = rms(VEL_ERR_E);
RMS_VEL_ERR_D = rms(VEL_ERR_D);


fprintf('SDR 속도 추정 결과');
fprintf('\n------------------------------------\n');
fprintf('[평균] 속도 오차(N)    : %f[m]\n',MEAN_VEL_ERR_N);
fprintf('[평균] 속도 오차(E)    : %f[m]\n',MEAN_VEL_ERR_E);
fprintf('[평균] 속도 오차(D)    : %f[m]\n',MEAN_VEL_ERR_D);
fprintf('[RMSE] 속도 오차(N)    : %f[m]\n',RMS_VEL_ERR_N);
fprintf('[RMSE] 속도 오차(E)    : %f[m]\n',RMS_VEL_ERR_E);
fprintf('[RMSE] 속도 오차(D)    : %f[m]\n',RMS_VEL_ERR_D);
fprintf('------------------------------------\n');


%% 속도 오차(Reference - SDR) Plot

figure;
subplot(3, 1, 1);
plot(TIME, VEL_ERR_N);
title('Velocity Error N');
xlabel('Time');
ylabel('Error [m/s]');
xlim([36, max(TIME)]); 

subplot(3, 1, 2);
plot(TIME, VEL_ERR_E);
title('Velocity Error E');
xlabel('Time');
ylabel('Error [m/s]');
xlim([36, max(TIME)]); 

subplot(3, 1, 3);
plot(TIME, VEL_ERR_D);
title('Velocity Error D');
xlabel('Time');
ylabel('Error [m/s]');
xlim([36, max(TIME)]); 


%% SDR 항법 결과(속도) Plot

VEL_N = SDR_VEL_TEMP(:, 2);
VEL_E = SDR_VEL_TEMP(:, 3);
VEL_D = SDR_VEL_TEMP(:, 4);

line_width = 2; 

% velocity - N
figure;
plot(TIME, VEL_N, 'linewidth', line_width);
grid on;
title('North Velocity');
xlabel('Time');
ylabel('V_N [m/s]');
ylim([-15, 15]);
xlim([36, max(TIME)]);

% velocity - E
figure;
plot(TIME, VEL_E, 'linewidth', line_width);
grid on;
title('East Velocity');
xlabel('Time');
ylabel('V_E [m/s]');
ylim([-15, 15]);
xlim([36, max(TIME)]);

% velocity - D
figure;
plot(TIME, VEL_D, 'linewidth', line_width);
grid on;
title('Down Velocity');
xlabel('Time');
ylabel('V_D [m/s]');
ylim([-15, 15]);
xlim([36, max(TIME)]);

