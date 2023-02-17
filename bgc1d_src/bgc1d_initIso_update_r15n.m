function bgc = bgc1d_initIso_update_r15n(bgc,tr)

 % % % % % % % % % % % % % % % % % % % % % % % % % %
 % This function updates isotope ratios 
 % % % % % % % % % % % % % % % % % % % % % % % % % %

 %%%%%  ratios  %%%%%
 %%%%%%%%%%%%%%%%%%%%

 %% Initialize isotopic ratios %%
 %bgc.r15norg = nan(size(bgc.zgrid));
 bgc.r15no3 = nan(size(bgc.nz));
 bgc.r15no2 = nan(size(bgc.nz));
 bgc.r15nh4 = nan(size(bgc.nz));
 bgc.r15n2o = nan(size(bgc.nz));
 bgc.r15n2oA = nan(size(bgc.nz));
 bgc.r15n2oB = nan(size(bgc.nz));

 %%%%%%%% Corg %%%%%%%
 %idx = tr.poc~=0;
 %ii = dNiso('d15N', bgc.d15norg, 'i14N', 1);
 %bgc.r15norg = ii.i15N./ii.N;     % 15N/(14N+15N) of organic nitrogen

 %%%%%%%% NO3 %%%%%%%%
 idx=(tr.no3+tr.i15no3~=0);
 bgc.r15no3(idx) = tr.i15no3(idx)./(tr.no3(idx)+tr.i15no3(idx));     % 15N/(14N+15N) of nitrate
 bgc.r14no3(idx) = (tr.no3(idx))./(tr.no3(idx)+tr.i15no3(idx));
 bgc.r15no3(~idx) = 0.0;
 bgc.r14no3(~idx)=0.0;

 %%%%%%%% NO2 %%%%%%%%
 idx=(tr.no2+tr.i15no2~=0);
 bgc.r15no2(idx) = tr.i15no2(idx)./(tr.no2(idx)+tr.i15no2(idx));     % 15N/(14N+15N) of nitrite
 bgc.r14no2(idx) = (tr.no2(idx))./(tr.no2(idx)+tr.i15no2(idx)); 
 bgc.r15no2(~idx) = 1.0;
 bgc.r14no2(~idx) = 1.0;

 %%%%%%%% NH4 %%%%%%%%
 idx=(tr.nh4+tr.i15nh4~=0);
 bgc.r15nh4(idx) = tr.i15nh4(idx)./(tr.nh4(idx)+tr.i15nh4(idx)); % 15N/(14N+15N) of ammonia
 bgc.r14nh4(idx) = (tr.nh4(idx))./(tr.nh4(idx)+tr.i15nh4(idx));
 bgc.r15nh4(~idx) = 1.0;
 bgc.r14nh4(~idx) = 1.0;

 %%%%%%% N2O %%%%%%%%
 tr.i15n2o = tr.i15n2oA + tr.i15n2oB;
 idx=(tr.i15n2o+tr.n2o~=0);
 bgc.r15n2oA(idx) = tr.i15n2oA(idx)./(tr.n2o(idx)+tr.i15n2o(idx)); % 15N/(14N+15N) of nitrous oxide Alpha
 bgc.r15n2oA(~idx)=1.0;
 bgc.r15n2oB(idx) = tr.i15n2oB(idx)./(tr.n2o(idx)+tr.i15n2o(idx)); % 15N/(14N+15N) of nitrous oxide Beta
 bgc.r15n2oB(~idx)=1.0;
 bgc.r15n2o(idx) = tr.i15n2o(idx)./(tr.n2o(idx)+tr.i15n2o(idx));     % 15N/(14N+15N) of nitrous oxide
 bgc.r14n2o(idx) = (tr.n2o(idx))./(tr.n2o(idx)+tr.i15n2o(idx));
 bgc.r15n2o(~idx) = 1.0;
 bgc.r14n2o(~idx) = 1.0;

%eps=10^-23;
%
% %%%%%%%% Corg %%%%%%%
% ii = dNiso('d15N', bgc.d15norg, 'N', tr.poc.*bgc.NCrem);
% bgc.r15norg = (ii.i15N+eps)./(ii.N+eps);     % 15N/(14N+15N) of organic nitrogen
%
% %%%%%%%% NO3 %%%%%%%%
% bgc.r15no3 = (tr.i15no3+eps)./(tr.no3+eps);     % 15N/(14N+15N) of nitrate
%
% %%%%%%%% NO2 %%%%%%%%
% bgc.r15no2 = (tr.i15no2+eps)./(tr.no2+eps);     % 15N/(14N+15N) of nitrite
%
% %%%%%%%% NH4 %%%%%%%%
% bgc.r15nh4 = (tr.i15nh4+eps)./(tr.nh4+eps);     % 15N/(14N+15N) of ammonia
%
% %%%%%%%% N2O %%%%%%%%
% bgc.r15n2oA = (tr.i15n2oA+eps)./(tr.n2o+eps); % 15N/(14N+15N) of nitrous oxide Alpha
% bgc.r15n2oB = (tr.i15n2oB+eps)./(tr.n2o+eps); % 15N/(14N+15N) of nitrous oxide Beta
% i15n2o = tr.i15n2oA + tr.i15n2oB;
% bgc.r15n2o = (i15n2o+eps)./(tr.n2o+eps);     % 15N/(14N+15N) of nitrous oxide
