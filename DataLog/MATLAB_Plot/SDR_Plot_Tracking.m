
clc;
clear;
close all;

for i = [10]
    FILE_NAME = sprintf('../DataLog_Trk_GPS%d.txt',i);
   
    Acq_Trc_Data_log = load (FILE_NAME);
    Acq_Trc_Data_log = Acq_Trc_Data_log(1:end,:);
    % 1 : 시간 에 따른 상관값 출력
    % 2 : 시간에 따른 DLL Discriminator 출력
    % 3 : 시간에 따른 FLL Discriminator 출력
    % 4 : 시간에 따른 DLL LoopFilter 출력
    % 5 : 시간에  따른 FLL LoopFilter 출력
    % 6 : 시간에 따른 PLL Discriminator 출력
    % 7 : 시간에 따른 PLL LoopFilter 출력
    % 8 : I_P
    % 9 : Q_P
    % 10 : Z_E
    % 11 : Z_L
    % 12 : Code Incriment
    % 13 : Carrier Doppler Frequency

   
    
    Z_std = std(Acq_Trc_Data_log(1000:end,1));
    Z_mean = mean(Acq_Trc_Data_log(1000:end,1));
    DLL_Discri_std = std(Acq_Trc_Data_log(1000:end,2));
    FLL_Discri_std = std(Acq_Trc_Data_log(2500:end,3));
    PLL_Discri_std = std(Acq_Trc_Data_log(end-1000:end,6))*180/pi;
    
%% E,P,L Correlator result

    figure();
    plot(Acq_Trc_Data_log(1:end,10),'r');       %E
    hold on;
    plot(Acq_Trc_Data_log(1:end,1),'b');        %P
    plot(Acq_Trc_Data_log(1:end,11),'g');       %L
    hold off;
    grid();
    title('상관값 Z(E,P,L)');    ylabel('상관값 Z');    xlabel('time [msec]');
    TITLE_NAME = sprintf('DataLog DLL FLL PRN : %d\n',i);    title(TITLE_NAME);
    legend('E','P','L');
    
    
%% DLL Discriminator 

    figure();
    plot(Acq_Trc_Data_log(:,2));                       
    hold on;
    plot(zeros(1,length(Acq_Trc_Data_log(:,3))),'g');
    hold off;
    TITLE_NAME = sprintf('DLL Discriminator Output(STD : %f[chip]).txt',DLL_Discri_std);   title(TITLE_NAME); 
    ylabel('DLL Discriminator[chip]');    xlabel('time [msec]');   grid();
    legend('DLL Discri[chip]','DLL Zero[Chip]');
    

%% FLL Discriminator plot

    figure();
    plot(Acq_Trc_Data_log(:,3));                        % 3 : FLL_Discri
    hold on;
    hold off;
    TITLE_NAME = sprintf('FLL Discriminator Output(STD : %f[Hz]).txt',FLL_Discri_std);   title(TITLE_NAME); 
    ylabel('FLL Discriminator[Hz]');    xlabel('time [msec]');   grid();
    legend('FLL Discri[Hz]','FLL Filtered Discri[Hz]','FLL Zero[Hz]');
    
    
%% PLL Discriminator plot   
    figure();
    plot(Acq_Trc_Data_log(:,6));                       
    TITLE_NAME = sprintf('PLL Discriminator Output(STD : %f[Hz]).txt',PLL_Discri_std);   title(TITLE_NAME); 
    ylabel('PLL Discriminator[Hz]');    xlabel('time [msec]');   grid();
    legend('PLL Discri[Hz]','PLL Filtered Discri[Hz]','PLL Zero[Hz]');
    
%% I/Q Diagram plot 
    figure()                                           
    plot(Acq_Trc_Data_log(end-1000:end,8), Acq_Trc_Data_log(end-1000:end,9), '.');
    grid();
    title('I-Q Plot');
    xlabel('P_I');
    ylabel('P_Q');
    IQ_mean = mean(sqrt((Acq_Trc_Data_log(end-1000:end,8)).^2 + (Acq_Trc_Data_log(end-1000:end,9)).^2));
    axis([-1*IQ_mean*1.2 IQ_mean*1.2 -1*IQ_mean*1.2 IQ_mean*1.2]);
   
    
    
    figure()
    subplot(2,1,1);
    plot(atan2(Acq_Trc_Data_log(end-1000:end,9),Acq_Trc_Data_log(end-1000:end,8))*180/(pi));
    title('atan2(Q,I) Output');
    ylabel('atan2(Q,I) Output[deg]');
    xlabel('time[msec]');
    axis([1 1000 -180 +180]);   grid();
    subplot(2,1,2);
    plot(atan(Acq_Trc_Data_log(end-1000:end,9)./Acq_Trc_Data_log(end-1000:end,8))*180/(pi));
    title('atan(Q/I) (PLL Discri)');
    ylabel('atan(Q/I) (PLL Discri)[deg]');
    xlabel('time[msec]');
    axis([1 1000 -90 +90]); grid();   
    
 
    
end