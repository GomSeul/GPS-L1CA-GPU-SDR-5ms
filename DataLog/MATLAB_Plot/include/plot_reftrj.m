function plot_reftrj(tsel, data, varargin)
    %   Plot reference trajectory
    %
    %   plot_reftrj(tsel, data, varargin)
    %
    %   INPUTS
    %       tsel = select time line index
    %               1: start from 0 second, 2: original
    %       data = position, velocity, attitude data with time index
    %       other arguments = same as plot function
    %
    %   OUTPUTS
    %       none
    %
    if tsel == 1
        t = data(:, 1) - data(1);
    else
        t = data(:, 1);
    end

    % Create a figure for Latitude and Longitude
    figure;
    subplot(2, 1, 1);
    plot(t, data(:, 2), varargin{:});
    grid on;
    ylabel('Latitude (deg)');
    xlabel('Time (s)');
    title('Latitude');

    subplot(2, 1, 2);
    plot(t, data(:, 3), varargin{:});
    grid on;
    ylabel('Longitude (deg)');
    xlabel('Time (s)');
    title('Longitude');

    % Create separate figures for other parameters
    param_names = {'Altitude (m)', 'V_N (m/s)', 'V_E (m/s)', 'V_D (m/s)', 'Roll (deg)', 'Pitch (deg)', 'Heading (deg)'};
    
    for param_idx = 4:10
        figure;
        plot(t, data(:, param_idx), varargin{:});
        grid on;
        ylabel(param_names{param_idx - 3});
        xlabel('Time (s)');
        title(param_names{param_idx - 3});
    end
end
%{
function plot_reftrj(tsel, data, varargin)
    %   Plot reference trajectory
    %
    %   plot_reftrj(tsel, data, varargin)
    %
    %   INPUTS
    %       tsel = select time line index
    %               1: start from 0 second, 2: original
    %       data = position, velocity, attitude data with time index
    %       other arguments = same as plot function
    %
    %   OUTPUTS
    %       none
    %
    if tsel == 1
        t = data(:, 1) - data(1);
    else
        t = data(:, 1);
    end

    % Create a figure for Latitude and Longitude
    figure;
    subplot(2, 1, 1);
    plot(t, data(:, 2), varargin{:});
    grid on;
    ylabel('Latitude (deg)');
    xlabel('Time (s)');

    subplot(2, 1, 2);
    plot(t, data(:, 3), varargin{:});
    grid on;
    ylabel('Longitude (deg)');
    xlabel('Time (s)');

    % Create separate figures for other parameters
    param_names = {'Altitude (m)', 'V_N (m/s)', 'V_E (m/s)', 'V_D (m/s)', 'Roll (deg)', 'Pitch (deg)', 'Heading (deg)'};
    
    for param_idx = 4:10
        figure;
        plot(t, data(:, param_idx), varargin{:});
        grid on;
        ylabel(param_names{param_idx - 3});
        xlabel('Time (s)');
    end
end
%}
%{
function plot_reftrj(tsel,data,varargin)
%   Plot reference trajectory
%
%   plot_reftrj(tsel,data,fig_num,...)
%
%   INPUTS
%       tsel = select time line index
%               1: start from 0 second, 2: original
%       data = position, velocity, attitude data with time index
%       other arguments = same as plot function
%
%   OUTPUTS
%       none
%
%
if tsel == 1
    t = data(:,1) - data(1);
    
else
    t = data(:,1);
end

figure;

%axis([14 300 70 110]);
subplot(3,3,1);
plot(t,data(:,2),varargin{:});
grid; 
axis tight;
ylabel('Latitude (deg)');


subplot(3,3,2);
plot(t,data(:,3),varargin{:});
grid; 
axis tight;
ylabel('Longitude (deg)');

%axis([-100 100 -100 100]);
%set(gca,'xtick',[-100:20:100]);
%set(gca,'ytick',[-100:20:100]);

%axis([0 200 126.45 128.5]);



subplot(3,3,3);
plot(t,data(:,4),varargin{:});
grid; 
axis tight;
ylabel('Altitude (m)');


%axis([0 160 70 110]);
subplot(3,3,4);
plot(t,data(:,5),varargin{:});
grid; 
axis tight;
ylabel('V_N (m/s)');
subplot(3,3,5);
plot(t,data(:,6),varargin{:});
grid; 
axis tight;
ylabel('V_E (m/s)');
subplot(3,3,6);
plot(t,data(:,7),varargin{:}); 
grid; 
axis tight;
% ylim([-7 -5]);
ylabel('V_D (m/s)');
subplot(3,3,7);
plot(t,data(:,8),varargin{:});
grid; 
axis tight;
% ylim([40 60]);
ylabel('Roll (deg)');
xlabel('Time (s)');
subplot(3,3,8);
plot(t,data(:,9),varargin{:});
grid; 
axis tight;
% ylim([4 6]);
ylabel('Pitch (deg)');
xlabel('Time (s)');
subplot(3,3,9);
plot(t,data(:,10),varargin{:});
grid; 
axis tight;
ylabel('Heading (deg)');
xlabel('Time (s)');
end
%}