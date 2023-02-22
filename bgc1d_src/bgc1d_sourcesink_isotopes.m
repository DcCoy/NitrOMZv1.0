function [sms diag] =  bgc1d_sourcesink_isotopes(bgc,tr)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Specifies the biogeochemical sources and sinks for NitrOMZ 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

 % calculate atom fractions
 if bgc.RunIsotopes
     % Update 15N/N ratios
     bgc = bgc1d_initIso_update_r15n(bgc,tr);
 end

% For safety, reduces zero and negative variables to small value
% Technically non-zero is mostly important for POC
tmpvar = fields(tr);
epsn = 1e-24;
for indf=1:length(tmpvar)
	tmp = max(epsn,tr.(tmpvar{indf}));
	tr.(tmpvar{indf}) = tmp;
end

tr.i15n2o = tr.i15n2oA+tr.i15n2oB;
% % % % % % % % % % % %
% % % % J-OXIC  % % % %
% % % % % % % % % % % %
% mm1 = Michaelis-Menton hyperbolic growth (var / var * k) where k is
% concentration of var where growth rate is half its maximum value
%----------------------------------------------------------------------
% (1) Oxic Respiration rate (C-units, mmolC/m3/s):
%----------------------------------------------------------------------
% Respiration rate based on POC, modified by oxygen concentration
RemOx = bgc.Krem .* mm1(tr.o2,bgc.KO2Rem) .* tr.poc;

%----------------------------------------------------------------------
% (2) Ammonium oxidation (molN-units, mmolN/m3/s):
%----------------------------------------------------------------------
% AO based on oxygen and NH4 concentration (loss of NH4)
Ammox = bgc.KAo .*  mm1(tr.o2,bgc.KO2Ao) .*  mm1_Iso(tr.nh4,tr.i15nh4,bgc.KNH4Ao) ;

%----------------------------------------------------------------------
% (3) Nitrite oxidation (molN-units, mmolN/m3/s):
%----------------------------------------------------------------------
% NO based on oxygen and NO2 concentration (loss of NO2)
Nitrox = bgc.KNo .*  mm1(tr.o2,bgc.KO2No) .* mm1_Iso(tr.no2,tr.i15no2,bgc.KNO2No);

%----------------------------------------------------------------------
% (4) N2O and NO2 production by ammox
%  Yields: nondimensional. Units of N, not N2O (molN-units, mmolN/m3/s): 
%----------------------------------------------------------------------
% Some AMMOX goes to N2O, most to NO2
Y = n2o_yield(tr.o2, bgc);
Jn2o_ao = Ammox.* Y.n2o;
Jno2_ao = Ammox.* Y.no2;

% % % % % % % % % % % %
% % %   J-ANOXIC  % % % 
% % % % % % % % % % % %

%----------------------------------------------------------------------
% (5) Denitrification (C-units, mmolC/m3/s)
%----------------------------------------------------------------------
RemDen1 = bgc.KDen1.* mm1_Iso(tr.no3,tr.i15no3,bgc.KNO3Den1) .* fexp(tr.o2,bgc.KO2Den1) .* tr.poc;
RemDen2 = bgc.KDen2 .* mm1_Iso(tr.no2,tr.i15no2,bgc.KNO2Den2) .* fexp(tr.o2,bgc.KO2Den2) .* tr.poc;
RemDen3 = bgc.KDen3 .* mm1_Iso(tr.n2o,(tr.i15n2oA+tr.i15n2oB),bgc.KN2ODen3) .* fexp(tr.o2,bgc.KO2Den3) .* tr.poc;
RemDen4 = bgc.KDen4 .* mm1_Iso(tr.no3,tr.i15no3,bgc.KNO3Den4) .* fexp(tr.o2,bgc.KO2Den4) .* tr.poc;

%----------------------------------------------------------------------
% (6) Anaerobic ammonium oxidation (molN2-units, mmolN2/m3/s):
% Note Anammox is in units of N2, i.e. 2 x mmolN/m3/s
%----------------------------------------------------------------------
Anammox = bgc.KAx .* mm1_Iso(tr.nh4,tr.i15nh4,bgc.KNH4Ax) .* mm1_Iso(tr.no2,tr.i15no2,bgc.KNO2Ax) .* fexp(tr.o2,bgc.KO2Ax);

%----------------------------------------------------------------------
% (8)  Calculate SMS for each tracer (mmol/m3/s)
%---------------------------------------------------------------------- 
sms.o2           =  -bgc.OCrem .* RemOx - 1.5.*Ammox - 0.5 .* Nitrox;
sms.no3          =  Nitrox - bgc.NCden1 .* RemDen1- bgc.NCden4 .* RemDen4; 
sms.poc          =  -(RemOx + RemDen1 + RemDen2 + RemDen3 + RemDen4);
sms.po4          =  bgc.PCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3 + RemDen4);
sms.nh4          =  bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3 + RemDen4) - (Jn2o_ao + Jno2_ao) - Anammox; % CLK: change 'Jnn2o_hx' to 'Jn2o_ao'
sms.no2          =  Jno2_ao + bgc.NCden1 .* RemDen1 - bgc.NCden2 .* RemDen2 - Anammox - Nitrox; 
sms.n2           =  bgc.NCden3 .* RemDen3 + Anammox; % (mmol N2/m3/s, units of N2, not N)
sms.kpoc         = -(RemOx + RemDen1 + RemDen2 + RemDen3)./tr.poc;
sms.n2oind.ammox = 0.5 .* Jn2o_ao; % (mmol N2O/m3/s, units of N2O, not N)
sms.n2oind.den2  = 0.5 .* bgc.NCden2 .* RemDen2; % (mmol N2O/m3/s, units of N2O, not N)
sms.n2oind.den3  = - bgc.NCden3 .* RemDen3; % (mmol N2O/m3/s, units of N2O, not N)
sms.n2oind.den4  = 0.5 .* bgc.NCden4 .* RemDen4;

% calculate binomial probabilities
[p1nh4, p2nh4, p3nh4, p4nh4] = binomial(bgc.r15nh4, bgc.r15nh4);
[p1no2, p2no2, p3no2, p4no2] = binomial(bgc.r15no2, bgc.r15no2);
[p1no3, p2no3, p3no3, p4no3] = binomial(bgc.r15no3, bgc.r15no3);

sms.n2o = (p4nh4 .* sms.n2oind.ammox ...
     + p4no2 .* sms.n2oind.den2 ...
     + sms.n2oind.den3 ...
     + p4no3 .* sms.n2oind.den4);

 if bgc.RunIsotopes
	 % Update 15N/N ratios
	 bgc = bgc1d_initIso_update_r15n(bgc,tr);
	 % Calculate sources and sinks for 15N tracers
	 sms.i15no3 = bgc.r15no2 .* bgc.alpha_nitrox .* Nitrox ...
	 	    - bgc.r15no3 .* bgc.alpha_den1 .* bgc.NCden1 .* RemDen1 ...
            - bgc.r15no3 .* bgc.alpha_den2 .* bgc.NCden4 .* RemDen4;
	 sms.i15no2 = bgc.r15nh4 .* (bgc.alpha_ammox_no2 .* Jno2_ao) ...
	            + bgc.r15no3 .* bgc.alpha_den1 .* bgc.NCden1 .* RemDen1 ...
	    	    - bgc.r15no2 .* bgc.alpha_den2 .* bgc.NCden2 .* RemDen2 ...
	            - bgc.r15no2 .* bgc.alpha_ax_no2 .* Anammox ...
	            - bgc.r15no2 .* bgc.alpha_nitrox .* Nitrox ;
	 sms.i15nh4 = bgc.r15norg.* bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3 + RemDen4) ...
	            - bgc.r15nh4 .* (bgc.alpha_ammox_no2 .* Jno2_ao) ...
	 	        - bgc.r15nh4 .* (bgc.alpha_ammox_n2o .* Jn2o_ao) ...
		        - bgc.r15nh4 .* bgc.alpha_ax_nh4 .* Anammox;
	 % N2O indivisual SMS 

     sms.i15n2oA = p2nh4 .* bgc.alpha_ammox_n2oA .* sms.n2oind.ammox ...
                 + p2no2 * bgc.alpha_den2A  .* sms.n2oind.den2 ...
                 + p2no3 * bgc.alpha_den2A  .* sms.n2oind.den4 ...
                 + bgc.r15n2oA .* bgc.alpha_den3_Alpha .* sms.n2oind.den3;
     sms.i15n2oB = p3nh4 .* bgc.alpha_ammox_n2oB .* sms.n2oind.ammox ...
                 + p3no2 * bgc.alpha_den2B  .* sms.n2oind.den2 ...
                 + p2no3 * bgc.alpha_den2B  .* sms.n2oind.den4 ...
                 + bgc.r15n2oB .* bgc.alpha_den3_Beta .* sms.n2oind.den3;
 end
%---------------------------------------------------------------------- 
% (9) Here adds diagnostics, to be handy when needed
%---------------------------------------------------------------------- 
diag.RemOx      = RemOx;                                       % mmolC/m3/s
diag.Ammox      = Ammox;                                       % mmolN/m3/s
diag.Nitrox     = Nitrox;                                      % mmolN/m3/s
diag.Anammox    = Anammox;                                     % mmolN2/m3/s
diag.RemDen1    = RemDen1;                                     % mmolC/m3/s
diag.RemDen2    = RemDen2;                                     % mmolC/m3/s
diag.RemDen3    = RemDen3;                                     % mmolC/m3/s
diag.RemDen     = RemDen1 + RemDen2 + RemDen3;                 % mmolC/m3/s
diag.Jno2_ao    = Jno2_ao;                                     % mmolN/m3/s
diag.Jn2o_ao    = Jn2o_ao;                                     % mmolN/m3/s
diag.Jn2o_prod  = sms.n2oind.ammox  + sms.n2oind.den2;         % mmolN2O/m3/s
diag.Jn2o_cons  = sms.n2oind.den3;                             % mmolN2O/m3/s
diag.Jno2_prod  = Jno2_ao + bgc.NCden1 .* RemDen1;             % mmolN/m3/s
diag.Jno2_cons  = - bgc.NCden2 .* RemDen2 - Anammox - Nitrox;  % mmolN/m3/s
diag.kpoc       = -(RemDen1 -RemDen2-RemDen3-RemOx) ./ tr.poc; % 1/s
%---------------------------------------------------------------------- 
