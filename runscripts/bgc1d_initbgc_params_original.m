function bgc = bgc1d_initbgc_params_original(bgc)

 %%%%%%%%% Stochiometry %%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Organic matter form: C_aH_bO_cN_dP
 % Anderson and Sarmiento 1994 stochiometry
 bgc.stoch_a = 106.0;
 bgc.stoch_b = 175.0;
 bgc.stoch_c = 42.0;
 bgc.stoch_d = 16.0;

%%%%%%%% Ammonification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.Krem = 0.08/(86400);               % 0.08    % Max. remineralization rate (1/s)
 bgc.KO2Rem  = 0.5;                     % 4       % Half sat. constant for respiration  (mmolO2/m3) - Martens-Habbena 2009

 %%%%%% Ammonium oxidationn %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Ammox: NH4 --> NO2
 bgc.KAo = 0.04556/(86400);              % 0.045    Max. Ammonium oxidation rate (1/s) - Bristow 2017
 bgc.KNH4Ao  = 0.0272;                  % 0.1       % Half sat. constant for nh4 (mmolN/m3) - Peng 2016
 bgc.KO2Ao = 0.333;                     % 0.333+-0.130  % Half sat. constant for Ammonium oxidation (mmolO2/m3) - Bristow 2017

 %%%%%%%% Nitrite oxidationn %%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Nitrox: NO2 --> NO3
 bgc.KNo =  0.255/(86400);             % 0.256    Max. Nitrite oxidation rate (1/s) - Bristow 2017
 bgc.KNO2No  = 0.0272;                          % Dont know (mmolN/m3)
 bgc.KO2No = 0.778;                     % Half sat. constant of NO2 for Nitrite oxidation (mmolO2/m3) - Bristow 2017

 %%%%%%% Denitrification %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Denitrif1: NO3 --> NO2
 bgc.KDen1 = 0.08/2/(86400);            % Max. denitrif1 rate (1/s)
 bgc.KO2Den1 = 1.0;                     % O2 poisoning constant for denitrif1 (mmolO2/m3)
 bgc.KNO3Den1 = 0.5;                    % Half sat. constant of NO3 for denitrif1 (mmolNO3/m3)

 % Denitrif2: NO2 --> N2O
 bgc.KDen2 = 0.08/6/(86400);            % Max. denitrif2 rate (1/s)
 bgc.KO2Den2 = 0.3;                    % O2 poisoning constant for denitrif2 (mmolO2/m3)
 bgc.KNO2Den2 = 0.5;                 % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3)

 % Denitrif3: N2O --> N2
 bgc.KDen3 = 0.08/3/(86400);           % Max. denitrif3 rate (1/s)
 bgc.KO2Den3 = 0.0292;                  % O2 poisoning constant for denitrif3 (mmolO2/m3)
 bgc.KN2ODen3 = 0.02;                   % Half sat. constant of N2O for denitrification3 (mmolNO3/m3)
 bgc.altKDen3 = 0.148/86400*0.6;       % alternate max N2O consumption rate, based on earlier model

 % Denitrif4: NO3 --> N2O
 bgc.KDen4 = 0.09/6/(86400); %0.08/6/(86400);            % Max. denitrif2 rate (1/s)
 bgc.KO2Den4 = 0.3;                    % Frey et al., 2020 (Figure 4c): 1/0.175 = 5.71
 bgc.KNO3Den4 = 0.5;                    % Half sat. constant of NO2 for denitrification2 (mmolNO3/m3)
 %%%%%%%%%%% Anammox %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.KAx = 0.02/86400;                % Max. Anaerobic Ammonium oxidation rate (1/s) - Bristow 2017
 bgc.KNH4Ax  = 0.0274;                  % Half sat. constant of NH4 for anammox (mmolNH4/m3)
 bgc.KNO2Ax  = 0.5;                  % Half sat. constant of NO2 for anammox (mmolNO2/m3)
 bgc.KO2Ax = 0.886;                               % 1.0     %

%%%%%%%%%% Hybrid N2O production %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.yHy     = 0.138 ;
bgc.KNO2Hy  = 0.1;
bgc.KO2Hy   = 0.19;          % O2 poisoning constant for hybrid N2O production

 %%%%%% N2O prod via ammox %%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Parameters for calculation of N2O yields during Ammox and
 % nitrifier-denitrification (see n2o_yield.m).
 
 % Choose paramterization:
 % 'Ji': Ji et al 2018
 bgc.n2o_yield = 'Ji';

 if strcmp(bgc.n2o_yield, 'Ji')
        % Ji et al 2018
        bgc.Ji_a = 0.2;
        bgc.Ji_b  = 0.08;
 end


