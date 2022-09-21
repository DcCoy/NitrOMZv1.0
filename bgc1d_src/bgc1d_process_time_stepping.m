function [dt_vec time_vec hist_time_vec hist_time_ind hist_time] = bgc1d_process_time_stepping(dt,endTimey,histTimey);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Get time-stepping variables 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if length(dt)~=length(endTimey)
	error('Wrong dt vector specification');
end

% number of different time steps to use
Ndt = length(dt);

% Builds vector of dts
dt_vec = [];
time_vec = 0;
for indi=1:Ndt
	tmp_ndt = floor((endTimey(indi)*86400*365-time_vec(end))/dt(indi));
	dt_vec = [dt_vec repmat(dt(indi),[1 tmp_ndt])]; 
	time_vec = cumsum(dt_vec);
end
hist_time_vec = [histTimey*86400*365:histTimey*86400*365:endTimey(end)*86400*365];
hist_time_ind = findin(hist_time_vec,time_vec);

% Updates history time vector (years)
hist_time = time_vec(hist_time_ind)/86400/365;


