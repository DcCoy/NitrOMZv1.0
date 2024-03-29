function bgc = bgc1d_initboundary(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Set top and bottom boundary conditions for specific regions 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
switch bgc.region
case 'ETNP'
   bgc.poc_flux_top = -6.0/86400*1.00; % Bry cond. for top POC Flux            (mmolC/m2/s)
   bgc.o2_top  = 203.65;               % Bry cond. for surface Oxygen          (mmolO2/m3)
   bgc.o2_bot  = 33.0;                 % Bry cond. for deep Oxygen             (mmolO2/m3)
   bgc.no3_top = 0.21;                 % Bry cond. for surface Nitrate         (mmolNO3/m3)
   bgc.no3_bot = 46.0;                 % Bry cond. for deep nitrate            (mmolNO3/m3)
   bgc.po4_top = 0.28;                 % Bry cond. for surface phosphate       (mmolPO4/m3)
   bgc.po4_bot =  3.38;                % Bry cond. for deep phosphate          (mmolPO4/m3)
   bgc.n2o_top = 6.0/1000;             % Bry cond. for surface nitrous oxide   (mmolN2O/m3)
   bgc.n2o_bot = 35.45/1000;           % Bry cond. for deep nitrous oxide      (mmolN2O/m3)
   bgc.nh4_top = 0;                    % Bry cond. for surface ammonia         (mmolNH4/m3)
   bgc.nh4_bot = 0;                    % Bry cond. for deep ammonia            (mmolNH4/m3)
   bgc.no2_top = 0;                    % Bry cond. for surface nitrite         (mmolNO2/m3)
   bgc.no2_bot = 0;                    % Bry cond. for deep nitrite            (mmolNO2/m3)
   bgc.n2_top = 2.0;                   % Bry cond. for surface N2 excess       (mmolN2/m3)
   bgc.n2_bot = 0.8131;                % Bry cond. for deep N2 excess          (mmolN2/m3)
case 'ETSP'
   bgc.poc_flux_top = -7.5/86400*0.8*1.85; % Bry cond. for top POC Flux            (mmolC/m2/s)
   bgc.o2_top  = 225;                      % Bry cond. for surface Oxygen          (mmolO2/m3)
   bgc.o2_bot  = 77;                       % Bry cond. for deep Oxygen             (mmolO2/m3)
   bgc.no3_top = 2.813;                    % Bry cond. for surface Nitrate         (mmolNO3/m3)
   bgc.no3_bot = 42.5;                     % Bry cond. for deep nitrate            (mmolNO3/m3)
   bgc.po4_top = 0.82;                     % Bry cond. for surface phosphate       (mmolPO4/m3)
   bgc.po4_bot =  3.055;                   % Bry cond. for deep phosphate          (mmolPO4/m3)
   bgc.n2o_top = 13.0/1000;                % Bry cond. for surface nitrous oxide   (mmolN2O/m3)
   bgc.n2o_bot = 35.0/1000;                % Bry cond. for deep nitrous oxide      (mmolN2O/m3)
   bgc.nh4_top = 0.4;                      % Bry cond. for surface ammonia         (mmolNH4/m3)
   bgc.nh4_bot = 10^-23;                   % Bry cond. for deep ammonia            (mmolNH4/m3)
   bgc.no2_top = 0.15;                     % Bry cond. for surface nitrite         (mmolNO2/m3)
   bgc.no2_bot = 10^-23;                   % Bry cond. for surface nitrite         (mmolNO2/m3)
   bgc.n2_top = 2.0;                       % Bry cond. for surface N2 excess       (mmolN2/m3)
   bgc.n2_bot = 6.0;                       % Bry cond. for deep N2 excess          (mmolN2/m3)
otherwise
   error(['Region ' bgc.region ' not found']);
   return
end
