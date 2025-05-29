function [df,time_diff] = diff_data(method,d1,d2)
%   Get difference between two data by comparing two time vectors
%
%	[df,time_diff] = diff_data(method,d1,d2)
%
%   INPUTS
%       method = the method for synchronization
%               1: union
%               2: intersection
%               3: uniform
%               4: nearest
%               5: exact
%       d1, d2 = the data to be processed
%               The first column must be time index in seconds.
%
%   OUTPUTS
%       df = time tagged difference between data
%       time_diff = time difference`
%
%

if method ~= 5
    [ds1, ds2] = get_synchronized_data(method,d1,d2);
else
    ds1 = d1;
    ds2 = d2;
end
df = [ds1(:,1),ds1(:,2:end) - ds2(:,2:end)];
time_diff = ds1(:,1) - ds2(:,1);