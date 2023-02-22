function bgc = bgc1d_initIso_update_r15n(bgc,t)

 % % % % % % % % % % % % % % % % % % % % % % % % % %
 % This function updates isotope ratios 
 % % % % % % % % % % % % % % % % % % % % % % % % % %

 %%%%%  ratios  %%%%%
 %%%%%%%%%%%%%%%%%%%%

 %% Initialize isotopic ratios %%
 bgc.r15norg = nan(size(bgc.zgrid));
 bgc.r15no3 = nan(size(bgc.zgrid));
 bgc.r15no2 = nan(size(bgc.zgrid));
 bgc.r15nh4 = nan(size(bgc.zgrid));
 bgc.r15n2o = nan(size(bgc.zgrid));
 bgc.r15n2oA = nan(size(bgc.zgrid));
 bgc.r15n2oB = nan(size(bgc.zgrid));

 %%%%%%%% Corg %%%%%%%
 idx = t.poc~=0;
 ii = dNiso('d15N', bgc.d15norg, 'i14N', 1);
 bgc.r15norg = ii.i15N./ii.N;     % 15N/(14N+15N) of organic nitrogen

 %%%%%%%% NO3 %%%%%%%%
 idx=(t.no3+t.i15no3~=0);
 bgc.r15no3(idx) = t.i15no3(idx)./(t.no3(idx)+t.i15no3(idx));     % 15N/(14N+15N) of nitrate
 bgc.r14no3(idx) = (t.no3(idx))./(t.no3(idx)+t.i15no3(idx));
 bgc.r15no3(~idx) = 0.0;
 bgc.r14no3(~idx)=0.0;

 %%%%%%%% NO2 %%%%%%%%
 idx=(t.no2+t.i15no2~=0);
 bgc.r15no2(idx) = t.i15no2(idx)./(t.no2(idx)+t.i15no2(idx));     % 15N/(14N+15N) of nitrite
 bgc.r14no2(idx) = (t.no2(idx))./(t.no2(idx)+t.i15no2(idx)); 
 bgc.r15no2(~idx) = 1.0;
 bgc.r14no2(~idx) = 1.0;

 %%%%%%%% NH4 %%%%%%%%
 idx=(t.nh4+t.i15nh4~=0);
 bgc.r15nh4(idx) = t.i15nh4(idx)./(t.nh4(idx)+t.i15nh4(idx)); % 15N/(14N+15N) of ammonia
 bgc.r14nh4(idx) = (t.nh4(idx))./(t.nh4(idx)+t.i15nh4(idx));
 bgc.r15nh4(~idx) = 1.0;
 bgc.r14nh4(~idx) = 1.0;

 %%%%%%% N2O %%%%%%%%
 t.i15n2o = t.i15n2oA + t.i15n2oB;
 idx=(t.i15n2o+t.n2o~=0);
 bgc.r15n2oA(idx) = t.i15n2oA(idx)./(t.n2o(idx)+t.i15n2o(idx)); % 15N/(14N+15N) of nitrous oxide Alpha
 bgc.r15n2oA(~idx)=1.0;
 bgc.r15n2oB(idx) = t.i15n2oB(idx)./(t.n2o(idx)+t.i15n2o(idx)); % 15N/(14N+15N) of nitrous oxide Beta
 bgc.r15n2oB(~idx)=1.0;
 bgc.r15n2o(idx) = t.i15n2o(idx)./(t.n2o(idx)+t.i15n2o(idx));     % 15N/(14N+15N) of nitrous oxide
 bgc.r14n2o(idx) = (t.n2o(idx))./(t.n2o(idx)+t.i15n2o(idx));
 bgc.r15n2o(~idx) = 1.0;
 bgc.r14n2o(~idx) = 1.0;

%eps=10^-23;
%
% %%%%%%%% Corg %%%%%%%
% ii = dNiso('d15N', bgc.d15norg, 'N', t.poc.*bgc.NCrem);
% bgc.r15norg = (ii.i15N+eps)./(ii.N+eps);     % 15N/(14N+15N) of organic nitrogen
%
% %%%%%%%% NO3 %%%%%%%%
% bgc.r15no3 = (t.i15no3+eps)./(t.no3+eps);     % 15N/(14N+15N) of nitrate
%
% %%%%%%%% NO2 %%%%%%%%
% bgc.r15no2 = (t.i15no2+eps)./(t.no2+eps);     % 15N/(14N+15N) of nitrite
%
% %%%%%%%% NH4 %%%%%%%%
% bgc.r15nh4 = (t.i15nh4+eps)./(t.nh4+eps);     % 15N/(14N+15N) of ammonia
%
% %%%%%%%% N2O %%%%%%%%
% bgc.r15n2oA = (t.i15n2oA+eps)./(t.n2o+eps); % 15N/(14N+15N) of nitrous oxide Alpha
% bgc.r15n2oB = (t.i15n2oB+eps)./(t.n2o+eps); % 15N/(14N+15N) of nitrous oxide Beta
% i15n2o = t.i15n2oA + t.i15n2oB;
% bgc.r15n2o = (i15n2o+eps)./(t.n2o+eps);     % 15N/(14N+15N) of nitrous oxide
