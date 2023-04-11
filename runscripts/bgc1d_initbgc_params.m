function bgc = bgc1d_initbgc_params(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Set stoichiometric ratios and bgc parameters 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%%% Stochiometry %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Organic matter form: C_aH_bO_cN_dP
% Stoichiometric ratios: C:H:O:N:P = a:b:c:d:1
% Anderson and Sarmiento 1994 stochiometry
bgc.stoch_a = 106.0;
bgc.stoch_b = 175.0;
bgc.stoch_c = 42.0;
bgc.stoch_d = 16.0;

%%%%%%% Ammonification %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.Krem = 0.08/86400 ;            % 0.08    % Max. specific remineralization rate (1/s)
bgc.KO2Rem  = 0.046571922641370;   % 4       % Half sat. constant for respiration  (mmolO2/m3) - Martens-Habbena 2009 %v5.4: 0.5;

%%%%%% Ammonium oxidationn %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ammox: NH4 --> NO2
bgc.KAo = 0.1170/86400;  % 0.045  % Max. Ammonium oxidation rate (mmolN/m3/s) - Bristow 2017 %v5.4: 0.04556/(86400);
bgc.KNH4Ao  = 0.130;     % 0.1    % Half sat. constant for nh4 (mmolN/m3) - Peng 2016 %v5.4: 0.0272;
bgc.KO2Ao = 0.333;       % 0.333+-0.130  % Half sat. constant for Ammonium oxidation (mmolO2/m3) - Bristow 2017

%%%%%%% Nitrite oxidation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nitrox: NO2 --> NO3
%bgc.KNo = 0.05/86400;   % 0.256 % Max. Nitrite oxidation rate (mmolN/m3/s) - Bristow 2017 %v5.4: 0.255/(86400);
bgc.KNo = 0.9/86400; %0.465/86400;   % ETNP 2018 rates;
bgc.KNO2No = 0.5;       % 1.0 % Dont know (mmolN/m3) %v5.4: 0.0272;
bgc.KO2No = 0.778;      % Half sat. constant of NO2 for Nitrite oxidation (mmolO2/m3) - Bristow 2017

%%%%%%% Denitrification %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denitrif1: NO3 --> NO2
%bgc.KDen1 = 0.0215/86400;  % Max. denitrif1 rate (1/s) %v5.4: 0.08/2/(86400);
bgc.KDen1 = 0.005/86400; %0.033/86400;  % ETNP 2018 rates;
bgc.KO2Den1 = 5.0;         % 5.0 O2 poisoning constant for denitrif1 (mmolO2/m3) %v5.4: 1.0;
bgc.KNO3Den1 = 0.4;        % Half sat. constant of NO3 for denitrif1 (mmolNO3/m3) %v5.4: 0.5;

% Denitrif2: NO2 --> N2O
bgc.KDen2 = 0.0005/86400;   %0.08/86400;% Max. denitrif2 rate (1/s) %v5.4: 0.08/6/(86400);
bgc.KO2Den2 = 2.3;         % O2 poisoning constant for denitrif2 (mmolO2/m3) %v5.4: 0.3;
bgc.KNO2Den2 = 0.05;       % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3) %v5.4: 0.5;

% Denitrif3: N2O --> N2
bgc.KDen3 = 0.0455/86400;   % Max. denitrif3 rate (1/s) %v5.4: 0.08/3/(86400); %v5.4: 0.08/3/(86400);
bgc.KO2Den3 = 0.11;         % O2 poisoning constant for denitrif3 (mmolO2/m3) %v5.4: 0.0292;
bgc.KN2ODen3 = 0.2;         % Half sat. constant of N2O for denitrification3 (mmolNO3/m3) %v5.4: 0.02;

% Denitrif4: NO3 --> N2O
bgc.KDen4 = 0.003/6/(86400);            %0.08/6/(86400); % Max. denitrif2 rate (1/s)
bgc.KO2Den4 = 0.3;                    % Frey et al., 2020 (Figure 4c): 1/0.175 = 5.71
bgc.KNO3Den4 = 0.5;                    % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3)
%%%%%%%%%% Anammox %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.KAx     = 0.3900/86400; % Max. Anaerobic Ammonium oxidation rate (mmolN/m3/s) - Bristow 2017 %v5.4: 0.02/86400;
bgc.KNH4Ax  = 0.23;         % Half sat. constant of NH4 for anammox (mmolNH4/m3) %v5.4: 0.0274;
bgc.KNO2Ax  = 0.1;          % Half sat. constant of NO2 for anammox (mmolNO2/m3) %v5.4: 0.5;
bgc.KO2Ax   = 0.7;          % O2 poisoning constant for anammox (mmolO2/m3) %v5.4: 0.886;

%%%%%%%%%% Hybrid N2O production %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.yHy     = 0.138 ;
bgc.KNO2Hy  = 0.1;
bgc.KO2Hy   = 0.19;          % O2 poisoning constant for hybrid N2O production

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


