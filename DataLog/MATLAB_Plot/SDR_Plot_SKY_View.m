%% SKY_PLOT


dataPath = '../DataLog_ELAZ_Out.txt';

data = importdata(dataPath);

% ������ ����
satelliteNumbers = data(:, 1);
elevationAngles = data(:, 2);
azimuthAngles = data(:, 3);

% SKY PLOT 
figure;
polarplot(0, 90, 'HandleVisibility', 'off'); 

hold on;

for i = 1:size(data, 1)
    az = azimuthAngles(i); % Azimuth ��
    el = 90 - elevationAngles(i); 
    
    % Azimuth�� Elevation�� SKYPLOT�� ǥ�� 
    h = polarplot(deg2rad(az), el, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'r', 'LineWidth', 2);
    h.MarkerFaceColor = 'none'; 
    
    % ���� ��ȣ ǥ�� 
    dx = 0.1; 
    text(deg2rad(az) + dx, el, num2str(satelliteNumbers(i)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
end

% SKYPLOT ��Ÿ�ϸ�
title('SKYPLOT of Visible Satellites');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
ax.RLim = [0, 90];
ax.RTick = [0, 30, 60, 90];
ax.RTickLabel = {'90��', '60��', '30��', '0��'};
ax.GridLineStyle = '--';

hold off;
