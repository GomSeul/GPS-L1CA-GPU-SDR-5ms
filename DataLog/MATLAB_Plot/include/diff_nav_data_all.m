function [df,time_diff] = diff_nav_data_all(method,d1,d2)
%   Get difference between two complete navigation data by comparing the
%   time vectors
%
%	[df,time_diff] = diff_nav_data_all(method,d1,d2)
%
%    INPUTS
%       method = the method for synchronization
%               1: union
%               2: intersection
%               3: uniform
%               4: nearest
%               5: exact
%       d1, d2 = the complete navigation data to be processed
%               The first column must be time index in seconds.
%
%    OUTPUTS
%       df = time tagged difference between two data
%       time_diff = time difference
%
%

[~, m1] = size(d1);
[~, m2] = size(d2);

if m1 < 10 || m2 < 10
    error('Insufficiency of number of vectors.');
end

% Reference 궤적
% 첫번쨰 인자 의 위치,속도,가속도 값을 구조체 in1 배열에 저장
t1 = d1(:,1);
in1(1).d = [t1,d1(:,2:4)];
in1(2).d = [t1,d1(:,5:7)];
in1(3).d = [t1,d1(:,8:10)];

% SDR 출력 
t2 = d2(:,1);
in2(1).d = [t2,d2(:,2:4)];
in2(2).d = [t2,d2(:,5:7)];
in2(3).d = [t2,d2(:,8:10)];

% 시간 범위 체크 (겹치는 부분이 없다면 오류 메시지 출력)
if t1(end) < t2(1) || t2(end) < t1(1)
    error('There dose not exist the common data period.');
end

parfor n = 1:3
    [df(n).d,df(n).time_diff] = diff_nav_data(method,n,in1(n).d,in2(n).d);
end

time_diff = df(1).time_diff;
df = [df(1).d,df(2).d(:,2:end),df(3).d(:,2:end)];
