

% Reference ��ǥ(�泲��)
%LAT= 36.3641948317362;
%LON = 127.345633430938;
%HEI = 93.7988147013063;

% Reference ��ǥ(��ϴ� A)
%LAT= 36.625013353;
%LON = 127.457856770;
%HEI = 110.303974416;

% Reference ��ǥ(��ϴ� B)
LAT= 36.625010955;
LON = 127.457860799;
HEI = 110.289323424;

shift = 0;

FILE_NAME = sprintf('../rover_cleaned.txt');
%FILE_NAME = sprintf('../DataLog_Navi_Out.txt');    %GPU 2ant A ���� �����
%FILE_NAME = sprintf('../GPU_ublox_AMP_true/DataLog_Navi_Out.txt');
DataLog_Navi = load (FILE_NAME);
DataLog_Navi = DataLog_Navi(1:end,:);

DataLog_Navi_Second     = DataLog_Navi(5:end,1+shift);      % receiver system time

DataLog_Navi_Pos_x      = DataLog_Navi(5:end,3+shift);      % position x
DataLog_Navi_Pos_y      = DataLog_Navi(5:end,4+shift);      % position y
DataLog_Navi_Pos_z      = DataLog_Navi(5:end,5+shift);      % position z

DataLog_Navi_lat        = DataLog_Navi(5:end,6+shift);      % position latitude
DataLog_Navi_lon        = DataLog_Navi(5:end,7+shift);      % position longitude
DataLog_Navi_hgt        = DataLog_Navi(5:end,8+shift);      % position high

DataLog_Navi_CLK_offset = DataLog_Navi(5:end,11+shift);   % clock offset


%% ���� ��ǥ ����

% *ones(length(DataLog_Navi_Pos_x),1)�� �ش� ũ�⸸ŭ�� �迭 ����
lat0   = (LAT)  *ones(length(DataLog_Navi_Pos_x),1);
lon0   = (LON)  *ones(length(DataLog_Navi_Pos_x),1);
h0      = (HEI)   *ones(length(DataLog_Navi_Pos_x),1);

% ECEF ��ǥ(�泲��)
%ecefx = (-3119380.78) *ones(length(DataLog_Navi_Pos_x),1);
%ecefy = (4088014.2) *ones(length(DataLog_Navi_Pos_y),1);
%ecefz = (3760865.44) *ones(length(DataLog_Navi_Pos_z),1);

% ECEF ��ǥ(��ϴ� A)
% ecefx = (-3116920.88) *ones(length(DataLog_Navi_Pos_x),1);
% ecefy = (4068242.70) *ones(length(DataLog_Navi_Pos_y),1);
% ecefz = (3784142.98) *ones(length(DataLog_Navi_Pos_z),1);

% ECEF ��ǥ(��ϴ� A)
ecefx = (-3116921.259) *ones(length(DataLog_Navi_Pos_x),1);
ecefy = (4068242.596) *ones(length(DataLog_Navi_Pos_y),1);
ecefz = (3784142.756) *ones(length(DataLog_Navi_Pos_z),1);


%% LLH ��ǥ PLOT

% Latitude/Longitude plot
figure();
plot(DataLog_Navi_lat, DataLog_Navi_lon, '.');
hold on;
plot(lat0, lon0, 'o');
TITLE_NAME = sprintf('Lat/Lon Plot');   title(TITLE_NAME); 
ylabel('Longitude[deg]');    xlabel('Latitude [deg]');
grid();

% Time/Latitude plot
figure();
subplot(2,1,1);
plot(DataLog_Navi_Second, DataLog_Navi_lat, '.');   
hold on;
plot(DataLog_Navi_Second,lat0, '-');
TITLE_NAME = sprintf('Latitude Plot');   
xlabel('Time [s]');   ylabel('Latitude [deg]');  title(TITLE_NAME); 
axis tight;  
grid();


subplot(2,1,2);
plot(DataLog_Navi_Second, DataLog_Navi_lon, '.');   
hold on;
plot(DataLog_Navi_Second,lon0, '-');
TITLE_NAME = sprintf('Longitude plot');   
xlabel('Time [s]');   ylabel('Longitude [deg]');  title(TITLE_NAME); 
axis tight;
grid();

% ���۷��� ��ǥ(lat0,lon0,h0) �������� ENU ��ǥ ���
% ���۷��� ��ǥ���� ������ SDR_Data_E/N/U�� ����
[SDR_Data_E SDR_Data_N SDR_Data_U] = ecef2enu(DataLog_Navi_Pos_x, ...
                                            DataLog_Navi_Pos_y, ...
                                            DataLog_Navi_Pos_z, ...
                                            lat0,lon0,h0,wgs84Ellipsoid);

%% RMSE ���
MEAN_POS_ERR_E = mean(SDR_Data_E);
MEAN_POS_ERR_N = mean(SDR_Data_N);
MEAN_POS_ERR_U = mean(SDR_Data_U);
MEAN_POS_ERR_EN = mean(sqrt(MEAN_POS_ERR_E^2 + MEAN_POS_ERR_N^2)); 
MEAN_POS_ERR_ENU = mean(sqrt(MEAN_POS_ERR_E^2 + MEAN_POS_ERR_N^2 + MEAN_POS_ERR_U^2)); 

MEAN_DIS_ERR_abs_U = mean(abs(SDR_Data_U));
MEAN_DIS_ERR_EN = mean(sqrt(((SDR_Data_E).^2 + (SDR_Data_N).^2)));
MEAN_DIS_ERR_ENU = mean(sqrt( ((SDR_Data_E).^2 + (SDR_Data_N).^2 + (SDR_Data_U).^2) ));

RMS_ERR_ENU = sqrt(mean((SDR_Data_E).^2 + (SDR_Data_N).^2+ (SDR_Data_U).^2));
RMS_ERR_EN = sqrt(mean(SDR_Data_E.^2 + SDR_Data_N.^2));
RMS_ERR_U = sqrt(mean(SDR_Data_U.^2));

TWO_RMS_EN = 2 * RMS_ERR_EN;
TWO_RMS_U = 2 * RMS_ERR_U;
TWO_RMS_ENU = 2 * RMS_ERR_ENU; 


fprintf('SDR ENU ������� ��� ����');

fprintf('\n------------------------------------\n');
fprintf('[���] �Ÿ�(����)      : %f[m]\n',MEAN_DIS_ERR_EN);
fprintf('[���] �Ÿ�(����)      : %f[m]\n',MEAN_DIS_ERR_abs_U);
fprintf('[���] �Ÿ�(����+����) : %f[m]\n',MEAN_DIS_ERR_ENU);

fprintf('[���] ��ġ ����(E)    : %f[m]\n',MEAN_POS_ERR_E);
fprintf('[���] ��ġ ����(N)    : %f[m]\n',MEAN_POS_ERR_N);
fprintf('[���] ��ġ ����(U)    : %f[m]\n',MEAN_POS_ERR_U);
fprintf('[���] ��ġ ����(EN)   : %f[m]\n',MEAN_POS_ERR_EN);
fprintf('[���] ��ġ ����(ENU)  : %f[m]\n',MEAN_POS_ERR_ENU);

fprintf('[RMSE] ����(ENU)      : %f[m]\n',RMS_ERR_EN);
fprintf('[RMSE] ����(U)        : %f[m]\n',RMS_ERR_U);
fprintf('[RMSE] ����+����(ENU) : %f[m]\n',RMS_ERR_ENU);
fprintf('------------------------------------\n');
                                        
%% ECEF(X/Y/Z) ��ǥ PLOT

% ECEF(3D) �÷�
figure();
plot3(DataLog_Navi_Pos_x,DataLog_Navi_Pos_y,DataLog_Navi_Pos_z,'o'); 
hold on;
plot3(ecefx,ecefy,ecefz,'o'); 
grid();
xlabel('x[m]');   ylabel('y[m]');   zlabel('z[m]');   title('ECEF(3D) Plot');  axis('equal');


% ECEF(X/Y/Z) �÷�
figure();
subplot(3,1,1);
plot(DataLog_Navi_Second,DataLog_Navi_Pos_x,'.-');   grid();
xlabel('Time [s]');   ylabel('x[m]');   title('ECEF(X/Y/Z) Plot');
axis tight;

subplot(3,1,2);
plot(DataLog_Navi_Second,DataLog_Navi_Pos_y,'.-');   grid();
xlabel('Time [s]');   ylabel('y[m]');  
axis tight;

subplot(3,1,3);
plot(DataLog_Navi_Second,DataLog_Navi_Pos_z,'.-');   grid();
xlabel('Time [s]');   ylabel('z[m]');  
axis tight;


%% ENU ��ǥ�� PLOT

% ENU(EN) �÷�
figure();
plot(SDR_Data_E, SDR_Data_N, 'o');  hold on;
plot(MEAN_POS_ERR_E, MEAN_POS_ERR_N, 'ro', 'MarkerSize', 6, 'LineWidth', 2.5); % ��� ���� ǥ��

%{
ezplot_fun = sprintf('(x-%f)^2+(y-%f)^2-%f',Diff_SDR_Pos_E, Diff_SDR_Pos_N, Std_SDR_Pos_EN^2);
h = ezplot(ezplot_fun,[Diff_SDR_Pos_E-Std_SDR_Pos_EN Diff_SDR_Pos_E+Std_SDR_Pos_EN Diff_SDR_Pos_N-Std_SDR_Pos_EN Diff_SDR_Pos_N+Std_SDR_Pos_EN]);
set(h,'LineWidth',2);  %# Sets the line width to 2
colormap([1 0 0]);
%}

hold off;
text_plot = sprintf('ENU(EN) Plot\nMean(E,N):(%f, %f)[m]',MEAN_POS_ERR_E,MEAN_POS_ERR_N);
xlabel('E[m]'); ylabel('N[m]'); title(text_plot);
axis('equal');   
hold on;   grid();
axis([-10 10 -10 10]);
set(gca,'xtick',[-10:2.5:10]);
set(gca,'ytick',[-10:2.5:10]);

% ENU(U) �÷�
figure();
plot(DataLog_Navi_Second, SDR_Data_U + HEI, '.-');
grid();
text_plot = sprintf('ENU(U) Plot\nMean(U):%f[m]', MEAN_POS_ERR_U);
xlabel('Time [s]');
ylabel('U Hight[m]');
title(text_plot);
axis tight;
ylim([80, 100]); % (�泲��)Y�� ������ 80���� 100���� ����, (��ϴ�)Y�� ������ 90���� 110���� ����
hold on;

% ENU(3D)�÷�
figure();
plot3(SDR_Data_E,SDR_Data_N,SDR_Data_U,'o'); 
grid();
xlabel('E[m]');   ylabel('N[m]');   zlabel('U[m]');   title('ENU(3D) Plot');  axis('equal');