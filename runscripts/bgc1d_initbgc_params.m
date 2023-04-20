function bgc = bgc1d_initbgc_params(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Set stoichiometric ratios and bgc parameters 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%%% Stochiometry %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organic matter form: C_aH_bO_cN_dP
% Stoichiometric ratios: C:H:O:N:P = a:b:c:d:1
% Anderson 1995 stochiometry
%bgc.stoch_a = 106.0;
%bgc.stoch_b = 175.0;
%bgc.stoch_c = 42.0;
%bgc.stoch_d = 16.0;
% Teng and Primeau 2014 stochiometry
bgc.stoch_a = 83;
bgc.stoch_b = 137;
bgc.stoch_c = 33;
bgc.stoch_d = 13;
% Tanioka 2022 stochiometry
%bgc.stoch_a = 126;
%bgc.stoch_b = 175.0;
%bgc.stoch_c = 42.0;
%bgc.stoch_d = 19;

%%%%%%% Ammonification %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.Krem = 0.0332/86400 ; % Max. specific remineralization rate (1/s) - Tiano (2014) offshore max rate, based on POCmax = 3.2 µM
%bgc.Krem = 0.35/86400 ;  % Max. specific remineralization rate (1/s) - Tiano (2014) station m1 (coastal) max rate, based on POCmax = 3.2 µM
%bgc.KO2Rem  = 0.018 ;       % Half sat. constant for respiration  (mmolO2/m3) - Tiano (2014) ETNP minimum value
bgc.KO2Rem  = 0.113 ;       % Half sat. constant for respiration  (mmolO2/m3) - Tiano (2014) ETNP maximum value

%%%%%% Ammonium oxidation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ammox: NH4 --> NO2
%bgc.KAo = 0.045/86400;  % Max. Ammonium oxidation rate (mmolN/m3/s) - Bristow 2016 %v5.4: 0.04556/(86400);
%bgc.KAo = 0.00217/86400; % Max. Ammonium oxidation rate (mmolN/m3/s) - Frey 2023 offshore
%bgc.KAo = 0.0059/86400;  % Max. Ammonium oxidation rate (mmolN/m3/s) - Frey  2023 coastal
bgc.KAo = 0.0904/86400; % Max. Ammonium oxidation rate (mmolN/m3/s) - Travis 2023
%bgc.KAo = 0.00221/86400; % Max. Ammonium oxidation rate (mmolN/m3/s) - CLK ETNP 2018 offshore;
%bgc.KAo = 0.00468/86400; % Max. Ammonium oxidation rate (mmolN/m3/s) - CLK ETNP 2018 coastal;

%bgc.KNH4Ao  = 0.134;    % Half sat. constant for nh4 (mmolN/m3) - Martens-Habbena 2009 %v5.4: 0.0272;
%bgc.KNH4Ao  = 0.0272;    % Half sat. constant for nh4 (mmolN/m3) - Peng 2016
bgc.KNH4Ao  = 0.106;    % Half sat. constant for nh4 (mmolN/m3) - Frey 2023 offshore
%bgc.KNH4Ao  = 0.048;    % Half sat. constant for nh4 (mmolN/m3) - Frey 2023 coastal
%bgc.KNH4Ao  = 0.81;    % Half sat. constant for nh4 (mmolN/m3) - N. viennensis Straka 2019

%bgc.KO2Ao = 0.333;       % 0.333+-0.130  % Half sat. constant for Ammonium oxidation (mmolO2/m3) - Bristow 2016
bgc.KO2Ao = 0.16;       % Half sat. constant for Ammonium oxidation (mmolO2/m3) - Frey 2023 offshore
%bgc.KO2Ao = 1.36;       % Half sat. constant for Ammonium oxidation (mmolO2/m3) - Frey 2023 coastal
%bgc.KO2Ao = 2.8;       % Half sat. constant for Ammonium oxidation (mmolO2/m3) - N. viennensis Straka 2019

%%%%%%% Nitrite oxidation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nitrox: NO2 --> NO3
%bgc.KNo = 0.256/86400; % Max. Nitrite oxidation rate (mmolN/m3/s) - Bristow 2016 %v5.4: 0.255/(86400);
%bgc.KNo = 0.03864/86400; % Max. Nitrite oxidation rate (mmolN/m3/s) - Sun 2021a
%bgc.KNo = 0.186/86400; % Max. Nitrite oxidation rate (mmolN/m3/s) - Babbin 2020 ETNP coastal
bgc.KNo = 0.113/86400; % Max. Nitrite oxidation rate (mmolN/m3/s) - Babbin 2020 ETNP offshore
%bgc.KNo = 0.0874/86400; % Max. Nitrite oxidation rate (mmolN/m3/s) - Travis 2023
%bgc.KNo = 0.465/86400;   % Max. Nitrite oxidation rate (mmolN/m3/s) - CLK ETNP 2018 rates;
%bgc.KNo = 0.01/86400; % Max. Nitrite oxidation rate (mmolN/m3/s) for anoxic nitrite oxidation

%bgc.KNO2No = 6.7;       % Half sat. constant for nitrite oxidation (mmolN/m3) - Sun 2021 offshore %v5.4: 0.0272;
bgc.KNO2No = 0.8;       % Half sat. constant for nitrite oxidation (mmolN/m3) - Sun 2021 coastal

%bgc.KO2No = 0.778;      % Half sat. constant of NO2 for Nitrite oxidation (mmolO2/m3) - Bristow 2017
bgc.KO2No = 0.01;      % Half sat. constant of NO2 for Nitrite oxidation (mmolO2/m3) - Sun 2021 (INSIGNIFICANT in sun 2021)

%%%%%%% Denitrification %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate [POC]max to convert rates to rate constants
% Denitrif1: NO3 --> NO2
%bgc.KDen1 = 0.0215/86400;  % Max. denitrif1 rate (1/s) %v5.4: 0.08/2/(86400);
%bgc.KDen1 = 0.00063/86400; % Max. denitrif1 rate (1/s) - Nicole light & rates paper coastal
bgc.KDen1 =  0.0024/86400; % Max. denitrif1 rate (1/s) - Nicole light & rates paper offshore
%bgc.KDen1 = 0.0011/86400; %0.003/86400; % Max. denitrif1 rate (1/s) - CLK ETNP coastal
%bgc.KDen1 = 0.0019/86400; %0.003/86400; % Max. denitrif1 rate (1/s) - CLK ETNP offshore
bgc.KO2Den1 = 5.0;         % O2 poisoning constant for denitrif1 (mmolO2/m3) %v5.4: 1.0;
bgc.KNO3Den1 = 0.4;        % Half sat. constant of NO3 for denitrif1 (mmolNO3/m3) %v5.4: 0.5;
%bgc.KDen1 =  0.01/86400; %0.0024/86400; % Max. denitrif1 rate (1/s) for anoxic nitrite oxidation

% Denitrif2: NO2 --> N2O
%bgc.KDen2 = 0.08/86400;% Max. denitrif2 rate (1/s) %v5.4: 0.08/6/(86400);
%bgc.KDen2 = 0.0000587/86400; % Max. denitrif2 rate (1/s) - CLK ETNP coastal
bgc.KDen2 =  0.0006/86400; %0.00000782/86400; % Max. denitrif2 rate (1/s) - CLK ETNP offshore
bgc.KO2Den2 = 0.9434;         % O2 poisoning constant for denitrif2 (mmolO2/m3) - CLK ETNP 2018 rates (1/1.06)
bgc.KNO2Den2 = 0.05;       % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3) %v5.4: 0.5;

% Denitrif3: N2O --> N2
%bgc.KDen3 = 0.0455/86400;   % Max. denitrif3 rate (1/s) %v5.4: 0.08/3/(86400);
%bgc.KDen3 = 0.000196/86400; % Max. denitrif3 rate (1/s)  - Sun 2021b coastal
bgc.KDen3 =  0.01/86400; %0.000362/86400;  % Max. denitrif3 rate (1/s)  - Sun 2021b offshore
%bgc.KDen3 = 0.000155/86400; % Max. denitrif3 rate (1/s)  - Babbin 2015
bgc.KO2Den3 = 0.11;         % O2 poisoning constant for denitrif3 (mmolO2/m3) %v5.4: 0.0292;
bgc.KN2ODen3 = 0.334;         % Half sat. constant of N2O for denitrification3 (mmolNO3/m3) %v5.4: 0.02;

% Denitrif4: NO3 --> N2O
%bgc.KDen4 = 0.08/6/(86400);    % Max. denitrif2 rate (1/s)
%bgc.KDen4 = 0.000375/86400;      % Max. denitrif2 rate (1/s) - CLK ETNP coastal
bgc.KDen4 =  0.0015/(86400); %0.000166/(86400); % Max. denitrif2 rate (1/s) - CLK ETNP offshore
bgc.KO2Den4 = 3.3333;           %  O2 poisoning constant for denitrif4 (mmolO2/m3) - CLK ETNP 2018 rates (1/0.3)
%bgc.KO2Den4 = 5.71;            %  O2 poisoning constant for denitrif4 (mmolO2/m3)- Frey et al., 2020 (Figure 4c)
bgc.KNO3Den4 = 0.5;             % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3)
%%%%%%%%%% Anammox %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bgc.KAx     = 0.0062/86400; % Max. Anaerobic Ammonium oxidation rate (mmolN/m3/s) - Bristow 2017 Bay of Bengal %v5.4: 0.02/86400;
%bgc.KAx     = 0.0262/86400; % Max. Anaerobic Ammonium oxidation rate (mmolN/m3/s) - Babbin 2020 ETNP coastal
bgc.KAx     = 0.0132/86400; % Max. Anaerobic Ammonium oxidation rate (mmolN/m3/s) - Babbin 2020 ETNP offshore
%bgc.KAx     = 0.0335/86400; % Max. Anaerobic Ammonium oxidation rate (mmolN/m3/s) - Karthäuser 2021 ETSP coastal
%bgc.KAx     = 0.0083/86400; % Max. Anaerobic Ammonium oxidation rate(mmolN/m3/s) - Karthäuser 2021 ETSP offshore
bgc.KNH4Ax  = 0.23; %3.0;         % Half sat. constant of NH4 for anammox (mmolNH4/m3) - Zhang 2020 %v5.4: 0.0274;
bgc.KNO2Ax  = 0.2; %0.1;          % Half sat. constant of NO2 for anammox (mmolNO2/m3) - Zhang 2020 %v5.4: 0.5;
bgc.KO2Ax   = 0.092;          % O2 poisoning constant for anammox (mmolO2/m3) - Straka 2019 %v5.4: 0.886;

%%%%%%%%%% Hybrid N2O production %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.yHy     = 0.20 ;          % maximum yield for hybrid N2O production
bgc.KNO2Hy  = 0.01;
bgc.KO2Hy   = 1.3514;          % O2 poisoning constant for hybrid N2O production

%%%%% N2O prod via ammox %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for calculation of N2O yields during Ammox and
% nitrifier-denitrification (see n2o_yield.m).

% Choose paramterization:
% 'Ji': Ji et al 2018
% 'Yang': Simon Yang 2019 optimized results
bgc.n2o_yield = 'Ji'; %'Yang';

switch bgc.n2o_yield
case 'Ji'
	% Ji et al 2018
	bgc.Ji_a = 0.2; 	% non-dimensional
	bgc.Ji_b  = 0.08; 	% 1/(mmolO2/m3)
case 'Yang'
	% S. Yang March 15 2019 optimized results 
	bgc.Ji_a = 0.30; 	% non-dimensional
	bgc.Ji_b  = 0.10; 	% 1/(mmolO2/m3)
otherwise
	error(['N2O yield ' bgc.n2o_yield  ' case not found']);
end


