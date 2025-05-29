function [df,time_diff] = diff_nav_data(method,sel,d1,d2)
%   Get difference between two navigation data by comparing the time
%   vectors
%
%	[df,time_diff] = diff_nav_data(method,sel,d1,d2)
%
%    INPUTS
%       method = the method for synchronization
%               1: union
%               2: intersection
%               3: uniform
%               4: nearest
%               5: exact
%       sel = navigation data type
%               1: position
%               2: velocity
%               3: attitude
%       d1, d2 = the navigation data to be processed
%               The first column must be time index in seconds.
%
%    OUTPUTS
%       df = time tagged difference between two data
%       time_diff = time difference
%
%

switch sel
    case 1
        if method ~= 5
            [ds1,ds2] = get_synchronized_data(method,d1,d2);
        else
            ds1 = d1;
            ds2 = d2;
        end
        df = convert_position_llh2ned(ds1(:,2:4),ds2(:,2:4));
        df = [ds1(:,1),df];
        time_diff = ds1(:,1) - ds2(:,1);
    case 2
        [df,time_diff] = diff_data(method,d1,d2);
    case 3
        % roll
        d1(:,2) = att_fix(d1(:,2));
        d2(:,2) = att_fix(d2(:,2));
        % heading
        d1(:,4) = att_fix(d1(:,4));
        d2(:,4) = att_fix(d2(:,4));
        [df,time_diff] = diff_data(method,d1,d2);
        % roll
        df(:,2) = att_fix(df(:,2));
        % heading
        df(:,4) = att_fix(df(:,4));
end
end

function y = att_fix(x)
% heading
y = x;
for k = 1:length(x)
    if y(k) > 180
        y(k) = y(k) - 360;
    elseif y(k) < -180
        y(k) = y(k) + 360;
    end
end
end