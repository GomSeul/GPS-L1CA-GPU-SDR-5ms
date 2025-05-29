
%% Acquisition plot (Peak Ratio)

    FILE_NAME = sprintf('../DataLog_Acq_Threshold.txt');
    Acq_Trc_Data_log = load (FILE_NAME);
    
figure();

for i = 1:32; 
    tY(i) = Acq_Trc_Data_log(i,1);
	tX(i) = i;   
end

bar(tX,tY); 

    title('Acquisition Result');	
	xlabel('PRN');
    ylabel('Peak Ratio');


