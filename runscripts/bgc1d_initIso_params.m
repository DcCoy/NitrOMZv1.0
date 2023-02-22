function bgc = bgc1d_initIso_params(bgc)

 %%%%%%%% General %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 bgc.IsoThreshold = 10^-10; % Threshold concentration for all N species exept n2o below which we do not 
 %			   compute deltas. For n2o bgc.IsoThreshold/1000 is used.

 %%%%%%%% Ammonification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 bgc.d15norg = 7.0;               % d15N of organic matter (permil)


 %%%%%% Ammonium oxidation %%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Epsilon for ...
 % .. ammox -- NH4-->NO2: 
 bgc.eps_ammox_no2 = 25.0;          % Enrichment factors during ammox (permil)
 % .. ammox -- NH4-->N2O:
 %bgc.eps_ammox_n2o = 47.0;          % this is more reflective of AOB - less predominant than AOA
 bgc.eps_ammox_n2oA = -21.3;        % Santoro et al., 2011 (adjusted based on δ15NH4+ value)
 bgc.eps_ammox_n2oB = 9.0;          % Santoro et al., 2011 (adjusted based on δ15NH4+ value)
 bgc.eps_ammox_n2o = (bgc.eps_ammox_n2oA + bgc.eps_ammox_n2oB)./ 2;
 % .. nitrif-denitrif -- NH4-->NO2:
 bgc.eps_nden_no2 = 15.0;           % Enrichment factors during nitrifier-denitrification (permil)
 % .. nitrif-denitrif -- NH4-->N2O:
 %bgc.eps_nden_n2o = 58.0;           % Frame and Casciotti [2010]
 bgc.eps_nden_n2oA = 63.5;           % Frame and Casciotti [2010]
 bgc.eps_nden_n2oB = 52.5;           % Frame and Casciotti [2010]
 bgc.eps_nden_n2o = (bgc.eps_nden_n2oA + bgc.eps_nden_n2oB) ./ 2;

 % SP for ...
 % .. ammox -- NH4-->N2O:
 bgc.n2oSP_ammox = 31.0;           % Site preference during N2O prod. via NH2OH (permil)
 % .. nitrif-denitrif -- NH4-->N2O:
 bgc.n2oSP_nden = -11.0;           % Site preference during nitrifier-denitrification (permil)


 %%%%%% Nitrite oxidation %%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Epsilon for nitrox -- NO2-->NO3
 bgc.eps_nitrox = -12.8;          % Enrichment factors during nitrox (permil)   

 %%%%%%% Denitrification %%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % epsilon for den1 -- NO3-->NO2
 bgc.eps_den1 = 15.0;          % Enrichment factors during denitrification 1 (permil)
 % epsilon for den2 -- NO2-->N2O
 %bgc.eps_den2 = 22.0;          % Enrichment factors during denitrification 2 (permil)
 bgc.eps_den2A = 0.0;          % Enrichment factors during denitrification 2 (permil)
 bgc.eps_den2B = 22.0;          % Enrichment factors during denitrification 2 (permil)
 bgc.eps_den2 = (bgc.eps_den2A + bgc.eps_den2B) ./ 2;
 % epsilon for den3 -- N2O-->N2
 bgc.eps_den3_Alpha = 20.0;          % Enrichment factor of N2OA during denitrification 3 (permil)
 bgc.eps_den3_Beta = 2.0;          % Enrichment factors of N2OB during denitrification 3 (permil)

 % SP for ..
 % .. den2 -- NO2-->N2O:
 bgc.n2oSP_den2 = -1.0;           % Site preference during denitrification 2 (permil)

 %%%%%%%%%%% Anammox %%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Epsilon for Anammox -- NO2-->N2
 bgc.eps_ax_no2 = 16.0;          % Enrichment factors during anammox (permil)
 % Epsilon for Anammox -- NH4-->N2
 bgc.eps_ax_nh4 = 25.0;          % Enrichment factors during anammox (permil)

 bgc = bgc1d_initIso_Dep_params(bgc);
