clc;
clear;
close all;

D2R = pi/180;
R2D = 180/pi;

%% REFERENCE ���� ������ ����

% ���Ͽ��� double������ �����͸� �о�ͼ� A ������ ���� + len�� 16���� ������ ������ ����
TRA_bin = fopen('2_Motion_scenario_80_stop_and_turn.bin');
A = fread(TRA_bin,'double');
len = length(A)/16;               % �������� ���� ��� 
ref_temp_1 = zeros(len,16);       %  ũ�Ⱑ (len, 16)�� 2���� �迭 ���� �� 0���� �ʱ�ȭ

% ù ��° ������ 0.1�� ��� + 40�̻��� ��쿡�� ����
for i = 1:len
    idx = (i-1)*16 + 1;
    if mod(A(idx),0.1) == 0  && A(idx)
        ref_temp_1(i,:) = A(idx:idx+15);
    end
end

% 0�� �� ����
ref_temp_1(~any(ref_temp_1,2),:) = [];
ref_temp_1 = round(ref_temp_1*10^7)/10^7;   % round to 7 decimal places
fclose(TRA_bin);

% ref_temp_1 �迭�� 1~7, 11~13���� �����Ͽ� REF_TRAJECTORY�� ����
REF_TRAJECTORY_1 = [ref_temp_1(:,1:7), ref_temp_1(:,11:13)];

%% ���۷��� ���� Plot
plot_3d_reftrj(REF_TRAJECTORY_1(:,2:4),1,'m','m');
plot_reftrj(1,REF_TRAJECTORY_1);


%% SDR �׹� ��� ����

% �ؽ�Ʈ ������ �Ǽ��� ���·� �迭�� ����
fileID = fopen('../Datalog_Navi_LLH_Out.txt','r');
sdr_temp = textscan(fileID,'%f %f %f %f','Delimiter','\t');
fclose(fileID);

sdr_temp = cell2mat(sdr_temp);
sdr_temp = round(sdr_temp*10^7)/10^7; % round to 7 decimal places


matching_rows = [];

matching_rows = [];
for row_ref = 1:size(ref_temp_1, 1)
    matching_indices = find(ismember(sdr_temp(:,1), ref_temp_1(row_ref, 1)));
    if ~isempty(matching_indices)
        matching_rows = [matching_rows, row_ref];
    end
end

ref_temp_2 = ref_temp_1(matching_rows, :);


% SDR �׹� ����� ������ �ð��� �����͸� ���۷��� ���� �迭�� ����
SDR_TRA = [sdr_temp(:,1:4),sdr_temp(:,1:4),sdr_temp(:,1:2)];
REF_TRAJECTORY_2 = [ref_temp_2(1:size(ref_temp_2, 1), 1:7), ref_temp_2(1:size(ref_temp_2, 1), 11:13)];


% ���� ���
nav_err = diff_nav_data_all(5,REF_TRAJECTORY_2,SDR_TRA);

%% �׹� ��� ���

nav_err_N = nav_err(:,2);
nav_err_E = nav_err(:,3);
nav_err_D = nav_err(:,4);

MEAN_POS_ERR_N = mean(nav_err_N);
MEAN_POS_ERR_E = mean(nav_err_E);
MEAN_POS_ERR_D = mean(nav_err_D);

MEAN_POS_ERR_NE = mean(sqrt(MEAN_POS_ERR_E^2 + MEAN_POS_ERR_N^2)); 
MEAN_POS_ERR_NED = mean(sqrt(MEAN_POS_ERR_E^2 + MEAN_POS_ERR_N^2 + MEAN_POS_ERR_D^2)); 

MEAN_DIS_ERR_abs_D = mean(abs(nav_err_D));
MEAN_DIS_ERR_NE = mean(sqrt(((nav_err_N).^2 + (nav_err_E).^2)));
MEAN_DIS_ERR_NED = mean(sqrt( ((nav_err_N).^2 + (nav_err_E).^2 + (nav_err_D).^2) ));

RMS_ERR_NED = sqrt(mean((nav_err_N).^2 + (nav_err_E).^2+ (nav_err_D).^2));
RMS_ERR_NE = sqrt(mean(nav_err_N.^2 + nav_err_E.^2));
RMS_ERR_D = sqrt(mean(nav_err_D.^2));

TWO_RMS_NE = 2 * RMS_ERR_NE;
TWO_RMS_U = 2 * RMS_ERR_D;
TWO_RMS_NED = 2 * RMS_ERR_NED; 



fprintf('SDR NED ������� ��� ����');

fprintf('\n------------------------------------\n');
fprintf('[���] �Ÿ�(����)      : %f[m]\n',MEAN_DIS_ERR_NE);
fprintf('[���] �Ÿ�(����)      : %f[m]\n',MEAN_DIS_ERR_abs_D);
fprintf('[���] �Ÿ�(����+����) : %f[m]\n',MEAN_DIS_ERR_NED);

fprintf('[���] ��ġ ����(E)    : %f[m]\n',MEAN_POS_ERR_E);
fprintf('[���] ��ġ ����(N)    : %f[m]\n',MEAN_POS_ERR_N);
fprintf('[���] ��ġ ����(D)    : %f[m]\n',MEAN_POS_ERR_D);
fprintf('[���] ��ġ ����(NE)   : %f[m]\n',MEAN_POS_ERR_NE);
fprintf('[���] ��ġ ����(NED)  : %f[m]\n',MEAN_POS_ERR_NED);

fprintf('[RMSE] ����(NED)      : %f[m]\n',RMS_ERR_NE);
fprintf('[RMSE] ����(D)        : %f[m]\n',RMS_ERR_D);
fprintf('[RMSE] ����+����(NED) : %f[m]\n',RMS_ERR_NED);
fprintf('------------------------------------\n');

%% NED�� �Ÿ� ���� Plot

figure;
plot(nav_err(:,1), nav_err(:,2:4));
xlabel('Time [sec]'); ylabel('Distance [m]');
legend('X', 'Y', 'Z');
title('Navigation Error');
grid on;
axis tight;



