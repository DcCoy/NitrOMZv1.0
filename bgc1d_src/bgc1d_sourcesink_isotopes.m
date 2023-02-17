function [sms diag] =  bgc1d_sourcesink_isotopes(bgc,tr)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Specifies the biogeochemical sources and sinks for NitrOMZ 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% handle zero isotopes
bgc = bgc1d_initIso_update_r15n(bgc,tr);

% For safety, reduces zero and negative variables to small value
% Technically non-zero is mostly important for POC
tmpvar = fields(tr);
epsn = 1e-24;
for indf=1:length(tmpvar) % try not adjusting isotopes to zero
	tmp = max(epsn,tr.(tmpvar{indf}));
	tr.(tmpvar{indf}) = tmp;
end

% Constant d15N of organic matter = 7 per mil
r15org = convert_delta(7, "d15N");

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
%Ammox = bgc.KAo .*  mm1(tr.o2,bgc.KO2Ao) .*  mm1(tr.nh4,bgc.KNH4Ao) ;
Ammox = bgc.KAo ./0.013 .* tr.nh4 .* mm1(tr.o2,bgc.KO2Ao);

%----------------------------------------------------------------------
% (3) Nitrite oxidation (molN-units, mmolN/m3/s):
%----------------------------------------------------------------------
% NO based on oxygen and NO2 concentration (loss of NO2)
%Nitrox = bgc.KNo .*  mm1(tr.o2,bgc.KO2No) .* mm1(tr.no2,bgc.KNO2No);
Nitrox = bgc.KNo ./3.2 * tr.no2 .*  mm1(tr.o2,bgc.KO2No);
Nitrox15 = bgc.KNo ./3.2 ./bgc.alpha_nitrox * tr.i15no2 .*  mm1(tr.o2,bgc.KO2No);

%----------------------------------------------------------------------
% (4) N2O and NO2 production by ammox
%  Yields: nondimensional. Units of N, not N2O (molN-units, mmolN/m3/s): 
%----------------------------------------------------------------------
% Some AMMOX goes to N2O, most to NO2
Y = n2o_yield(tr.o2, bgc);
Jn2o_ao = bgc.KAo ./0.013 .* tr.nh4 .* mm1(tr.o2,bgc.KO2Ao) .* Y.n2o;
Jno2_ao = bgc.KAo ./0.013 .* tr.nh4 .* mm1(tr.o2,bgc.KO2Ao) .* Y.no2;
Jn2o_ao15 = max(0, bgc.KAo ./0.013 ./bgc.alpha_ammox_n2o .* tr.i15nh4 .* mm1(tr.o2,bgc.KO2Ao) .* Y.n2o);
Jno2_ao15 = max(0, bgc.KAo ./0.013 ./bgc.alpha_ammox_no2 .* tr.i15nh4 .* mm1(tr.o2,bgc.KO2Ao) .* Y.no2);
    
% % % % % % % % % % % %
% % %   J-ANOXIC  % % % 
% % % % % % % % % % % %

%----------------------------------------------------------------------
% (5) Denitrification (C-units, mmolC/m3/s)
%----------------------------------------------------------------------
%RemDen1 = bgc.KDen1 .* mm1(tr.no3,bgc.KNO3Den1) .* fexp(tr.o2,bgc.KO2Den1) .* tr.poc;
RemDen1 = bgc.KDen1 ./48 .* tr.no3 .* fexp(tr.o2,bgc.KO2Den1) .* tr.poc;
RemDen1_15 = max(0, bgc.KDen1 ./48 ./bgc.alpha_den1 .* tr.i15no3 .* fexp(tr.o2,bgc.KO2Den1) .* tr.poc);

%RemDen2 = bgc.KDen2 .* mm1(tr.no2,bgc.KNO2Den2) .* fexp(tr.o2,bgc.KO2Den2) .* tr.poc;
RemDen2 = bgc.KDen2 ./3.2 .* tr.no2 .* fexp(tr.o2,bgc.KO2Den2) .* tr.poc;
RemDen2_15 = bgc.KDen2 ./3.2 ./bgc.alpha_den2 .* tr.i15no2 .* fexp(tr.o2,bgc.KO2Den2) .* tr.poc;

%RemDen3 = bgc.KDen3 .* mm1(tr.n2o,bgc.KN2ODen3) .* fexp(tr.o2,bgc.KO2Den3) .* tr.poc;
RemDen3 = bgc.KDen3 ./0.14 .* tr.n2o .* fexp(tr.o2,bgc.KO2Den3) .* tr.poc;

%----------------------------------------------------------------------
% (6) Anaerobic ammonium oxidation (molN2-units, mmolN2/m3/s):
% Note Anammox is in units of N2, i.e. 2 x mmolN/m3/s
%----------------------------------------------------------------------
%Anammox = bgc.KAx .* mm1(tr.nh4,bgc.KNH4Ax) .* mm1(tr.no2,bgc.KNO2Ax) .* fexp(tr.o2,bgc.KO2Ax);
%Anammox = max(0, bgc.KAx ./3.2 ./0.14 .* tr.nh4 .* tr.no2 .* fexp(tr.o2,bgc.KO2Ax));
%Anammox15nh4 = max(0, bgc.KAx ./3.2 ./0.14 ./bgc.alpha_ax_nh4 .* tr.i15nh4 .* tr.i15no2 .* fexp(tr.o2,bgc.KO2Ax));
%Anammox15no2 = max(0, bgc.KAx ./3.2 ./0.14 ./bgc.alpha_ax_no2 .* tr.i15nh4 .* tr.i15no2 .* fexp(tr.o2,bgc.KO2Ax));

Anammox = bgc.KAx ./0.14 .* tr.nh4 .* fexp(tr.o2,bgc.KO2Ax); % NH4-only anammox term for 15NH4-only troubleshooting
Anammox15nh4 = bgc.KAx ./0.14 ./bgc.alpha_ax_nh4 .* tr.i15nh4 .* fexp(tr.o2,bgc.KO2Ax);
Anammox15no2 = bgc.KAx ./0.14 ./bgc.alpha_ax_no2 .* tr.i15nh4 .* fexp(tr.o2,bgc.KO2Ax);

%----------------------------------------------------------------------
% (8)  Calculate SMS for each tracer (mmol/m3/s)
%---------------------------------------------------------------------- 
sms.o2           =  -bgc.OCrem .* RemOx - 1.5.*Ammox - 0.5 .* Nitrox;
sms.no3          =  Nitrox - bgc.NCden1 .* RemDen1; 
sms.poc          =  -(RemOx + RemDen1 + RemDen2 + RemDen3);
sms.po4          =  bgc.PCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3);
sms.nh4          =  bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3) - (Jn2o_ao + Jno2_ao) - Anammox; % CLK: change 'Jnn2o_hx' to 'Jn2o_ao'
sms.no2          =  Jno2_ao + bgc.NCden1 .* RemDen1 - bgc.NCden2 .* RemDen2 - Anammox - Nitrox; 
sms.n2           =  bgc.NCden3 .* RemDen3 + Anammox; % (mmol N2/m3/s, units of N2, not N)
sms.kpoc         = -(RemOx + RemDen1 + RemDen2 + RemDen3)./tr.poc;
sms.n2oind.ammox = 0.5 .* Jn2o_ao; % (mmol N2O/m3/s, units of N2O, not N)
sms.n2oind.den2  = 0.5 .* bgc.NCden2 .* RemDen2; % (mmol N2O/m3/s, units of N2O, not N)
sms.n2oind.den3  = - bgc.NCden3 .* RemDen3; % (mmol N2O/m3/s, units of N2O, not N)
sms.n2o          = (sms.n2oind.ammox + sms.n2oind.den2 + sms.n2oind.den3);

if bgc.RunIsotopes
sms.i15no3 =  RemOx * 0;
sms.i15no2 =  r15org .* bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3)- (Jn2o_ao15 + Jno2_ao15) - Anammox15no2;
sms.i15nh4 =  r15org .* bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3)- (Jn2o_ao15 + Jno2_ao15) - Anammox15nh4;
sms.i15n2oA = RemOx * 0;
sms.i15n2oB = RemOx * 0;
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
