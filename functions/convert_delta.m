function R = convert_delta(d, i)

R15std = 0.00367647; % air N2
R18std = 0.00200517; % VSMOW

if i == "d15N"
    R = ((d/1000)+1)*R15std;
elseif i == "d18O"
    R = ((d/1000)+1)*R18std;
else
    disp('Please enter i as d15N or d18O');
    R = -999;    
end