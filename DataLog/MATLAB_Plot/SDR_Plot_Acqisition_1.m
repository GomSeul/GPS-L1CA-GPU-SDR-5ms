
    
%% Acquisition plot (Full Search)

f_s = 25000000;
sample_size_1ms = (f_s/1000);

range_dopp_freq = -11*10^3:500:11*10^3; % (Hz) 도플러 주파수 범위
range_shift_code = [0:1:sample_size_1ms-1]/sample_size_1ms*1023; % (chip) 코드 천이 범위


 for i = [2]  
     
    FILE_NAME = sprintf('../Datalog_Acq_Full_FFT_GPS%d.txt',i);
    Z_acq = load (FILE_NAME);

    [X,Y] = meshgrid(range_shift_code, range_dopp_freq);
    figure()
    mesh(X,Y,Z_acq);
    grid();
  
    axis([0 1023 -11*10^3 11*10^3 0 max(max(Z_acq))*1.1]);
    title_data = sprintf('GPS L1 C/A : %d',i);
    title(title_data);
    xlabel('code shift (0.5 chip)');
    ylabel('doppler freq (500Hz)');
    zlabel('Z mag');
    
end
