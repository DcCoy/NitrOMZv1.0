 function [sms diag] =  bgc1d_sourcesink_isos(bgc,t); 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % iNitrOMZ v1.0 - Simon Yang  - April 2019
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Specifies the biogeochemical sources and sinks for OMZ nutrient/POC/N2O model
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Make sure we have no negative concentration
 if bgc.RunIsotopes
     % Update 15N/N ratios
     bgc = bgc1d_initIso_update_r15n(bgc,t);
 end
 ff=fields(t);
 for f = 1 : length(ff)
 	t.(ff{f})(t.(ff{f})<0)=0;
 end

 % % % % % % % % % % % %
 % % % % J-OXIC  % % % %
 % % % % % % % % % % % %
     %----------------------------------------------------------------------
     % (1) Oxic Respiration rate (C-units):
     %----------------------------------------------------------------------
     RemOx = bgc.Krem .* mm1(t.o2,bgc.KO2Rem) .* t.poc;

 if ~bgc.RunIsotopes 
     %----------------------------------------------------------------------
     % (2) Ammonium oxidation (molN-units):
     %----------------------------------------------------------------------
     Ammox = bgc.KAo .*  mm1(t.o2,bgc.KO2Ao) .*  mm1(t.nh4,bgc.KNH4Ao) ;

     %----------------------------------------------------------------------
     % (3) Nitrite oxidation (molN-units):
     %----------------------------------------------------------------------
     Nitrox = bgc.KNo .*  mm1(t.o2,bgc.KO2No) .* mm1(t.no2,bgc.KNO2No);

     %----------------------------------------------------------------------
     % (4) N2O and NO2 production by ammox and nitrifier-denitrif (molN-units): 
     %----------------------------------------------------------------------
     Y = n2o_yield_nden(t.o2, bgc);
     % via NH2OH
     Jnn2o_hx = Ammox .* Y.nn2o_hx_nh4;
     Jno2_hx = Ammox.* Y.no2_hx_nh4; %Ammox .* (Y.no2_hx_nh4+Y.nn2o_nden_nh4);
     % via NH4->NO2->N2O
     Jnn2o_nden = Ammox .* Y.nn2o_nden_nh4; %0;
     
     % % % % % % % % % % % %
     % % %   J-ANOXIC  % % % 
     % % % % % % % % % % % %

     %----------------------------------------------------------------------
     % (5) Denitrification (C-units)
     %----------------------------------------------------------------------
     RemDen1 = bgc.KDen1.* mm1(t.no3,bgc.KNO3Den1) .* fexp(t.o2,bgc.KO2Den1) .* t.poc;
     RemDen2 = bgc.KDen2 .* mm1(t.no2,bgc.KNO2Den2) .* fexp(t.o2,bgc.KO2Den2) .* t.poc;
     RemDen3 = bgc.KDen3 .* mm1(t.n2o,bgc.KN2ODen3) .* fexp(t.o2,bgc.KO2Den3) .* t.poc;
     RemDen4 = bgc.KDen4 .* mm1(t.no3,bgc.KNO3Den4) .* fexp(t.o2,bgc.KO2Den4) .* t.poc;

     %----------------------------------------------------------------------
     % (6) Anaerobic ammonium oxidation (molN-units):
     %----------------------------------------------------------------------
     Anammox = bgc.KAx .* mm1(t.nh4,bgc.KNH4Ax) .* mm1(t.no2,bgc.KNO2Ax) .* fexp(t.o2,bgc.KO2Ax);
	
     % Non-isotope case
     bgc.r14no3 = 1.0;
     bgc.r14no2 = 1.0;
     bgc.r14nh4 = 1.0;
     bgc.r14n2o = 1.0;
 else
     t.i15n2o = t.i15n2oA+t.i15n2oB;
     %----------------------------------------------------------------------
     % (2) Ammonium oxidation (molN-units):
     %----------------------------------------------------------------------
     Ammox = bgc.KAo .*  mm1(t.o2,bgc.KO2Ao) .*  mm1_Iso(t.nh4,t.i15nh4,bgc.KNH4Ao) ;

     %----------------------------------------------------------------------
     % (3) Nitrite oxidation (molN-units):
     %----------------------------------------------------------------------
     Nitrox = bgc.KNo .*  mm1(t.o2,bgc.KO2No) .* mm1_Iso(t.no2,t.i15no2,bgc.KNO2No);

     %----------------------------------------------------------------------
     % (4) N2O and NO2 production by ammox and nitrifier-denitrif (molN-units):
     %----------------------------------------------------------------------
     Y = n2o_yield_nden(t.o2, bgc);
     % via NH2OH
     Jnn2o_hx =Ammox .* Y.nn2o_hx_nh4;
     Jno2_hx = Ammox.* Y.no2_hx_nh4;
     % via NH4->NO2->N2O
     Jnn2o_nden = Ammox .* Y.nn2o_nden_nh4;

     % % % % % % % % % % % %
     % % %   J-ANOXIC  % % %
     % % % % % % % % % % % %

     %----------------------------------------------------------------------
     % (5) Denitrification (C-units)
     %----------------------------------------------------------------------
     RemDen1 = bgc.KDen1.* mm1_Iso(t.no3,t.i15no3,bgc.KNO3Den1) .* fexp(t.o2,bgc.KO2Den1) .* t.poc;
     RemDen2 = bgc.KDen2 .* mm1_Iso(t.no2,t.i15no2,bgc.KNO2Den2) .* fexp(t.o2,bgc.KO2Den2) .* t.poc;
     RemDen3 = bgc.KDen3 .* mm1_Iso(t.n2o,t.i15n2o,bgc.KN2ODen3) .* fexp(t.o2,bgc.KO2Den3) .* t.poc;
     RemDen4 = bgc.KDen4 .* mm1_Iso(t.no3,t.i15no3,bgc.KNO3Den4) .* fexp(t.o2,bgc.KO2Den4) .* t.poc;

     %----------------------------------------------------------------------
     % (6) Anaerobic ammonium oxidation (molN-units):
     %----------------------------------------------------------------------
     Anammox = bgc.KAx .* mm1_Iso(t.nh4,t.i15nh4,bgc.KNH4Ax) .* mm1_Iso(t.no2,t.i15no2,bgc.KNO2Ax) .* fexp(t.o2,bgc.KO2Ax);
 end

 KRemOx = RemOx./t.poc;
 KRemDen1 = RemDen1./t.poc;
 KRemDen2 = RemDen2./t.poc;
 KRemDen3 = RemDen3./t.poc;

 %----------------------------------------------------------------------
 % (8)  Calculate SMS for each tracer
 %---------------------------------------------------------------------- 
 sms.o2   =  (-bgc.OCrem .* RemOx - 1.5.*Ammox - 0.5 .* Nitrox);
 sms.no3  =  (Nitrox - bgc.NCden1 .* RemDen1- bgc.NCden4 .* RemDen4);% .*bgc.r14no3;
 sms.poc  =  (- RemOx - (RemDen1+RemDen2+RemDen3));
 sms.po4  =  (+bgc.PCrem .* (RemOx) + bgc.PCden1 .* RemDen1 + bgc.PCden2 .* RemDen2 + bgc.PCden3 .* RemDen3+ bgc.PCden3 .* RemDen4);
 sms.nh4  =  (+ bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3 + RemDen4) - (Jnn2o_hx+Jno2_hx+Jnn2o_nden) - Anammox);% .*bgc.r14nh4;
 sms.no2  =  (Jno2_hx + bgc.NCden1 .* RemDen1 - bgc.NCden2 .* RemDen2 - Anammox - Nitrox);% .*bgc.r14no2;
 sms.n2   =  (bgc.NCden3 .* RemDen3 + Anammox);
 sms.kpoc = -(KRemOx + KRemDen1 + KRemDen2 + KRemDen3);
 % N2O individual SMSs
 sms.n2oind.ammox = 0.5 .* Jnn2o_hx;
 sms.n2oind.nden  = 0.5 .* Jnn2o_nden;
 sms.n2oind.den2  = 0.5 .* bgc.NCden2 .* RemDen2;
 sms.n2oind.den3  = - bgc.NCden3 .* RemDen3;
 sms.n2oind.den4 = 0.5 .* bgc.NCden4 .* RemDen4;
 %sms.n2oind.den3  = - bgc.altKDen3.* fexp(t.o2,bgc.KO2Den3) .*t.n2o; % with first-order rate law

 % N2O total SMS
 if ~bgc.RunIsotopes
 sms.n2o = (sms.n2oind.ammox ...
     + sms.n2oind.nden ...
     + sms.n2oind.den2 ...
     + sms.n2oind.den3 ...
     + sms.n2oind.den4);
 
 else
 % calculate binomial probabilities
 [p1nh4, p2nh4, p3nh4, p4nh4] = binomial(bgc.r15nh4, bgc.r15nh4);
 [p1no2, p2no2, p3no2, p4no2] = binomial(bgc.r15no2, bgc.r15no2);
 [p1no3, p2no3, p3no3, p4no3] = binomial(bgc.r15no3, bgc.r15no3);

 % N2O total SMS
 sms.n2o = (p4nh4 .* sms.n2oind.ammox ...
     + p4nh4 .* sms.n2oind.nden ...
     + p4no2 .* sms.n2oind.den2 ...
     + sms.n2oind.den3 ...
     + p4no3 .* sms.n2oind.den4);

 
	 % Update 15N/N ratios
	 bgc = bgc1d_initIso_update_r15n(bgc,t);
     %disp(bgc.r15nh4)
	 % Calculate sources and sinks for 15N tracers
	 sms.i15no3 = bgc.r15no2 .* bgc.alpha_nitrox .* Nitrox ...
	 	    - bgc.r15no3 .* bgc.alpha_den1 .* bgc.NCden1 .* RemDen1 ...
            - bgc.r15no3 .* bgc.alpha_den2 .* bgc.NCden4 .* RemDen4;
	 sms.i15no2 = bgc.r15nh4 .* (bgc.alpha_ammox_no2 .* Jno2_hx) ...
	            + bgc.r15no3 .* bgc.alpha_den1 .* bgc.NCden1 .* RemDen1 ...
	    	    - bgc.r15no2 .* bgc.alpha_den2 .* bgc.NCden2 .* RemDen2 ...
	            - bgc.r15no2 .* bgc.alpha_ax_no2 .* Anammox ...
	            - bgc.r15no2 .* bgc.alpha_nitrox .* Nitrox ;
	 sms.i15nh4 = bgc.r15norg.* bgc.NCrem .* (RemOx + RemDen1 + RemDen2 + RemDen3 + RemDen4) ...
	            - bgc.r15nh4 .* (bgc.alpha_ammox_no2 .* Jno2_hx) ...
	 	        - bgc.r15nh4 .* (bgc.alpha_ammox_n2o .* Jnn2o_hx + bgc.alpha_nden_n2o .* Jnn2o_nden) ...
		        - bgc.r15nh4 .* bgc.alpha_ax_nh4 .* Anammox;
	 % N2O indivisual SMS 

     sms.i15n2oA = p2nh4 .* bgc.alpha_ammox_n2oA .* sms.n2oind.ammox ...
                 + p2nh4 .* bgc.alpha_nden_n2oA .* sms.n2oind.nden ...
                 + p2no2 * bgc.alpha_den2A  .* sms.n2oind.den2 ...
                 + p2no3 * bgc.alpha_den2A  .* sms.n2oind.den4 ...
                 + bgc.r15n2oA .* bgc.alpha_den3_Alpha .* sms.n2oind.den3;
     sms.i15n2oB = p3nh4 .* bgc.alpha_ammox_n2oB .* sms.n2oind.ammox ...
                 + p3nh4 .* bgc.alpha_nden_n2oB .* sms.n2oind.nden ...
                 + p3no2 * bgc.alpha_den2B  .* sms.n2oind.den2 ...
                 + p3no3 * bgc.alpha_den2B  .* sms.n2oind.den4 ...
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
diag.RemDen4    = RemDen4;                                     % mmolC/m3/s
diag.RemDen     = RemDen1 + RemDen2 + RemDen3 + RemDen4;       % mmolC/m3/s
diag.Jno2_ao    = Jno2_hx;                                     % mmolN/m3/s
diag.Jn2o_ao    = Jnn2o_hx + Jnn2o_nden;                       % mmolN/m3/s
diag.Jn2o_prod  = sms.n2oind.ammox + sms.n2oind.nden + sms.n2oind.den2 + sms.n2oind.den4;         % mmolN2O/m3/s
diag.Jn2o_cons  = sms.n2oind.den3;                             % mmolN2O/m3/s
diag.Jno2_prod  = Jno2_hx + bgc.NCden1 .* RemDen1;             % mmolN/m3/s
diag.Jno2_cons  = - bgc.NCden2 .* RemDen2 - Anammox - Nitrox;  % mmolN/m3/s
diag.kpoc       = -(RemDen1 -RemDen2-RemDen3-RemOx) ./ t.poc; % 1/s
%---------------------------------------------------------------------- 
if bgc.RunIsotopes
	sms.r15no3=bgc.r15no3;
	sms.r15no2=bgc.r15no2;
	sms.r15nh4=bgc.r15nh4;
	sms.r15norg=bgc.r15nh4;
	sms.r15n2o=bgc.r15n2o;
end