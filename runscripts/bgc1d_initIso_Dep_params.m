 function bgc = bgc1d_initIso_Dep_params(bgc)

 % % % % % % % % % % % % % % % % % % % % % % % % % % 
 % This function Calculates fractionation factors
 % from enrichment factors as well as 15N 
 % concentrations for boundary conditions
 % % % % % % % % % % % % % % % % % % % % % % % % % %

 %%%%%% Ammonium oxidation %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for ammox -- NH4-->NO2:
 bgc.alpha_ammox_no2 = 1.0 - bgc.eps_ammox_no2./1000.0;  % Fractionation factor during ammox (permil)
 % Alpha for ammox -- NH4-->N2O:
 bgc.alpha_ammox_n2o = 1.0 - bgc.eps_ammox_n2o./1000.0;  % Fractionation factor during N2O prod. via NH2OH (permil)
 bgc.alpha_ammox_n2oA = 1.0 - bgc.eps_ammox_n2oA./1000.0;  % Fractionation factor during N2O prod. via NH2OH (permil)
 bgc.alpha_ammox_n2oB = 1.0 - bgc.eps_ammox_n2oB./1000.0;  % Fractionation factor during N2O prod. via NH2OH (permil)
 % Alpha for nitrif-denitrif -- NH4-->NO2:
 bgc.alpha_nden_no2 = 1.0 - bgc.eps_nden_no2./1000.0;    % Fractionation factor during nitrifier-denitrification (permil)
 % Alpha for nitrif-denitrif -- NH4-->N2O:
 bgc.alpha_nden_n2o = 1.0 - bgc.eps_nden_n2o./1000.0;    % Fractionation factor during nitrifier-denitrification (permil)
 bgc.alpha_nden_n2oA = 1.0 - bgc.eps_nden_n2o./1000.0;    % Fractionation factor during nitrifier-denitrification (permil)
 bgc.alpha_nden_n2oB = 1.0 - bgc.eps_nden_n2o./1000.0;    % Fractionation factor during nitrifier-denitrification (permil)

 %%%%%% Nitrite oxidation %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for nitrox -- NO2-->NO3
 bgc.alpha_nitrox = 1.0 - bgc.eps_nitrox./1000.0;        % Fractionation factor during nitrox (permil)

 %%%%%%% Denitrification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for den1 -- NO3-->NO2
 bgc.alpha_den1 = 1.0 - bgc.eps_den1./1000.0;          % Fractionation factor during denitrification 1 (permil)
 % Alpha for den2 -- NO2-->N2O
 bgc.alpha_den2 = 1.0 - bgc.eps_den2./1000.0;          % Fractionation factor during denitrification 2 (permil)
 bgc.alpha_den2A = 1.0 - bgc.eps_den2A./1000.0;          % Fractionation factor during denitrification 2 (permil)
 bgc.alpha_den2B = 1.0 - bgc.eps_den2B./1000.0;          % Fractionation factor during denitrification 2 (permil)
 % Alpha for den3 -- N2OA-->N2
 bgc.alpha_den3_Alpha = 1.0 - bgc.eps_den3_Alpha./1000.0;          % Fractionation factor during denitrification 3 (permil)
  % Alpha for den3 -- N2OB-->N2
 bgc.alpha_den3_Beta = 1.0 - bgc.eps_den3_Beta./1000.0;          % Fractionation factor during denitrification 3 (permil)

 %%%%%%%%%%% Anammox %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for Anammox -- NO2-->N2
 bgc.alpha_ax_no2 = 1.0 - bgc.eps_ax_no2./1000.0;        % Fractionation factor during anammox (permil)
 % Alpha for Anammox -- NH4-->N2
 bgc.alpha_ax_nh4 = 1.0 - bgc.eps_ax_nh4./1000.0;        % Fractionation factor during anammox (permil)

 %%%%%%  bound. conditions %%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%% NO3 %%%%%%%%
 % Top boundary
 ii = dNiso('d15N', bgc.d15no3_top, 'N', bgc.no3_top);
 bgc.i15no3_top = ii.i15N;     % Concentration of 15N nitrate
 % Deep boundary
 ii = dNiso('d15N', bgc.d15no3_bot, 'N', bgc.no3_bot);
 bgc.i15no3_bot = ii.i15N;     % Concentraion of 15N nitrate

 %%%%%%%% NO2 %%%%%%%%
 % Top boundary
 ii = dNiso('d15N', bgc.d15no2_top, 'N', bgc.no2_top);
 bgc.i15no2_top = ii.i15N;     % Concentration of 15N nitrite
 % Deep boundary
 ii = dNiso('d15N', bgc.d15no2_bot, 'N', bgc.no2_bot);
 bgc.i15no2_bot = ii.i15N;     % Concentration of 15N nitrite

 %%%%%%%% NH4 %%%%%%%%
 % Top boundary
 ii = dNiso('d15N', bgc.d15nh4_top, 'N', bgc.nh4_top);
 bgc.i15nh4_top = ii.i15N;     % Concentration of 15N ammonia
 % Deep boundary
 ii = dNiso('d15N', bgc.d15nh4_bot, 'N', bgc.nh4_bot);
 bgc.i15nh4_bot = ii.i15N;     % Concentration of 15N ammonia

 %%%%%%%% N2OA %%%%%%%%
 % Top boundary
 ii = dNiso('d15N', bgc.d15n2oA_top, 'N', bgc.n2o_top);
 bgc.i15n2oA_top = ii.i15N;     % Concentration of 15N nitrous oxide Alpha
 % Deep boundary
 ii = dNiso('d15N', bgc.d15n2oA_bot, 'N', bgc.n2o_bot);
 bgc.i15n2oA_bot = ii.i15N;     % Concentration of 15N nitrous oxide Alpha

 %%%%%%%% N2OB %%%%%%%%
 % Top boundary
 ii = dNiso('d15N', bgc.d15n2oB_top, 'N', bgc.n2o_top);
 bgc.i15n2oB_top = ii.i15N;     % Concentration of 15N nitrous oxide Beta
 % Deep boundary
 ii = dNiso('d15N', bgc.d15n2oB_bot, 'N', bgc.n2o_bot);
 bgc.i15n2oB_bot = ii.i15N;     % Concentration of 15N nitrous oxide Beta

