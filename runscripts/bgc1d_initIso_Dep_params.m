 function bgc = bgc1d_initIso_Dep_params(bgc)

 % % % % % % % % % % % % % % % % % % % % % % % % % % 
 % This function Calculates fractionation factors
 % from enrichment factors as well as 15N 
 % concentrations for boundary conditions
 % % % % % % % % % % % % % % % % % % % % % % % % % %

 %%%%%% Ammonium oxidation %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CLK: changed alphas to formulation I'm more familiar with
 % alpha = eps/1000 + 1
 % Alpha for ammox -- NH4-->NO2:
 bgc.alpha_ammox_no2 = bgc.eps_ammox_no2./1000.0 + 1.0;  % Fractionation factor during ammox (permil)
 % Alpha for ammox -- NH4-->N2O:
 bgc.alpha_ammox_n2o = bgc.eps_ammox_n2o./1000.0 + 1.0;  % Fractionation factor during N2O prod. via NH2OH (permil)
 % Alpha for nitrif-denitrif -- NH4-->NO2:
 bgc.alpha_nden_no2 = bgc.eps_nden_no2./1000.0 + 1.0;    % Fractionation factor during nitrifier-denitrification (permil)
 % Alpha for nitrif-denitrif -- NH4-->N2O:
 bgc.alpha_nden_n2o = bgc.eps_nden_n2o./1000.0 + 1.0;    % Fractionation factor during nitrifier-denitrification (permil)

 %%%%%% Nitrite oxidation %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for nitrox -- NO2-->NO3
 bgc.alpha_nitrox = bgc.eps_nitrox./1000.0 + 1.0;        % Fractionation factor during nitrox (permil)

 %%%%%%% Denitrification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for den1 -- NO3-->NO2
 bgc.alpha_den1 = bgc.eps_den1./1000.0 + 1.0;          % Fractionation factor during denitrification 1 (permil)
 % Alpha for den2 -- NO2-->N2O
 bgc.alpha_den2 = bgc.eps_den2./1000.0 + 1.0;          % Fractionation factor during denitrification 2 (permil)
 % Alpha for den3 -- N2OA-->N2
 bgc.alpha_den3_Alpha = bgc.eps_den3_Alpha./1000.0 + 1.0;          % Fractionation factor during denitrification 3 (permil)
  % Alpha for den3 -- N2OB-->N2
 bgc.alpha_den3_Beta = bgc.eps_den3_Beta./1000.0 + 1.0;          % Fractionation factor during denitrification 3 (permil)

 %%%%%%%%%%% Anammox %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Alpha for Anammox -- NO2-->N2
 bgc.alpha_ax_no2 = bgc.eps_ax_no2./1000.0 + 1.0;        % Fractionation factor during anammox (permil)
 % Alpha for Anammox -- NH4-->N2
 bgc.alpha_ax_nh4 = bgc.eps_ax_nh4./1000.0 + 1.0;        % Fractionation factor during anammox (permil)

 %%%%%%  bound. conditions %%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %%%%%%%% NO3 %%%%%%%%
 % Top boundary
 %ii = dNiso('d15N', bgc.d15no3_top, 'N', bgc.no3_top);
 %bgc.i15no3_top = ii.i15N;
 % CLK: replaced dNiso with a simpler function
 bgc.i15no3_top = convert_delta(bgc.d15no3_top,'d15N')*bgc.no3_top;     % Concentration of 15N nitrate
 % Deep boundary
 bgc.i15no3_bot = convert_delta(bgc.d15no3_bot,'d15N')*bgc.no3_bot;     % Concentration of 15N nitrate

 %%%%%%%% NO2 %%%%%%%%
 % Top boundary
 bgc.i15no2_top = convert_delta(bgc.d15no2_top,'d15N')*bgc.no2_top;     % Concentration of 15N nitrite
 % Deep boundary
 bgc.i15no2_bot = convert_delta(bgc.d15no2_bot,'d15N')*bgc.no2_bot;     % Concentration of 15N nitrite

 %%%%%%%% NH4 %%%%%%%%
 % Top boundary
 bgc.i15nh4_top = convert_delta(bgc.d15nh4_top,'d15N')*bgc.nh4_top;     % Concentration of 15N ammonia
 % Deep boundary
 bgc.i15nh4_bot = convert_delta(bgc.d15nh4_bot,'d15N')*bgc.nh4_bot;     % Concentration of 15N ammonia

 %%%%%%%% N2OA %%%%%%%%
 % Top boundary
 bgc.i15n2oA_top = convert_delta(bgc.d15n2oA_top,'d15N')*bgc.n2o_top;     % Concentration of 15N nitrous oxide Alpha
 % Deep boundary
 bgc.i15n2oA_bot = convert_delta(bgc.d15n2oA_bot,'d15N')*bgc.n2o_bot;     % Concentration of 15N nitrous oxide Alpha

 %%%%%%%% N2OB %%%%%%%%
 % Top boundary
 bgc.i15n2oB_top = convert_delta(bgc.d15n2oB_top,'d15N')*bgc.n2o_top;     % Concentration of 15N nitrous oxide Beta
 % Deep boundary
 bgc.i15n2oB_bot = convert_delta(bgc.d15n2oB_bot,'d15N')*bgc.n2o_bot;     % Concentration of 15N nitrous oxide Beta
