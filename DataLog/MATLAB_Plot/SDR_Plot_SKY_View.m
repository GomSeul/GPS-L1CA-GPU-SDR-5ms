%% SKY_PLOT


dataPath = '../DataLog_ELAZ_Out.txt';

data = importdata(dataPath);

% 데이터 추출
satelliteNumbers = data(:, 1);
elevationAngles = data(:, 2);
azimuthAngles = data(:, 3);

% SKY PLOT 
figure;
polarplot(0, 90, 'HandleVisibility', 'off'); 

hold on;

for i = 1:size(data, 1)
    az = azimuthAngles(i); % Azimuth 값
    el = 90 - elevationAngles(i); 
    
    % Azimuth와 Elevation을 SKYPLOT에 표시 
    h = polarplot(deg2rad(az), el, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    h.MarkerFaceColor = 'none'; 
    
    % 위성 번호 표시 
    dx = 0.1; 
    text(deg2rad(az) + dx, el, num2str(satelliteNumbers(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end

% SKYPLOT 스타일링
title('SKYPLOT of Visible Satellites');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
ax.RLim = [0, 90];
ax.RTick = [0, 30, 60, 90];
ax.RTickLabel = {'90°', '60°', '30°', '0°'};
ax.GridLineStyle = '--';

hold off;
