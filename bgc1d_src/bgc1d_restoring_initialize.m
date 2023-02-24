function  bgc = bgc_initialize_restoring(bgc);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Interpolate profile by creating a cubic spline 
% NOTE: modify substituting [0;tconc;0] if flat-slope endings are needed
% WARNING: avoid extrapolation by using profiles from 0-4000 m (even if sparse)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

load([bgc.root,bgc.farfield_profiles]);
load([bgc.root,bgc.Tau_profiles]);

% Get farfield file
bgc.restore = [];
if  strcmp(bgc.region, 'ETNP')
   Farfield = farfield_ETNP_gridded;
elseif strcmp(bgc.region, 'ETSP')
   Farfield = farfield_ETSP_gridded;
end

% tauZvar
if bgc.tauZvar == 1
   rest.tauh = ((rest.currentZ_etnp./(bgc.Lh) + rest.kappaZ_etnp.*2./((bgc.Lh)^2)).^-1);
   tdepth = -abs(rest.depth);
   tconc  = rest.tauh;     
   ibad = find(isnan(tconc));
   tdepth(ibad) = [];
   tconc(ibad) = [];     
   bgc.restore.tauh_pre = spline(tdepth,tconc);
end

% PO4
if ~(bgc.PO4rest==0)
   tdepth = -abs(Farfield.zgrid);
   tconc  = Farfield.po4;
   ibad = find(isnan(tconc));
   tdepth(ibad) = [];
   tconc(ibad) = [];
   bgc.restore.po4_cs = spline(tdepth,tconc);
   bgc.restore.po4_cout = ppval(bgc.restore.po4_cs,bgc.zgrid);    
end

% NO3
if ~(bgc.NO3rest==0)
   tdepth = -abs(Farfield.zgrid);
   tconc  = Farfield.no3;
   ibad = find(isnan(tconc));
   tdepth(ibad) = [];
   tconc(ibad) = [];
   bgc.restore.no3_cs = spline(tdepth,tconc);
   bgc.restore.no3_cout = ppval(bgc.restore.no3_cs,bgc.zgrid);    
end

% O2
if ~(bgc.O2rest==0)
   tdepth = -abs(Farfield.zgrid);
   tconc  = Farfield.o2;
   ibad = find(isnan(tconc));
   tdepth(ibad) = [];
   tconc(ibad) = [];
   bgc.restore.o2_cs = spline(tdepth,tconc);
   bgc.restore.o2_cout = ppval(bgc.restore.o2_cs,bgc.zgrid);    
end

% N2O
if ~(bgc.N2Orest==0)
   tdepth = -abs(Farfield.zgrid);
   tconc  = Farfield.n2o;
   ibad = find(isnan(tconc));
   tdepth(ibad) = [];
   tconc(ibad) = [];
   bgc.restore.n2o_cs = spline(tdepth,tconc);
   bgc.restore.n2o_cout = ppval(bgc.restore.n2o_cs,bgc.zgrid);    
end

% NO2
if ~(bgc.NO2rest==0)
   tdepth = -abs(Farfield.zgrid);
   tconc  = Farfield.no2;
   ibad = find(isnan(tconc));
   tdepth(ibad) = [];
   tconc(ibad) = [];
   bgc.restore.no2_cs = spline(tdepth,tconc);
   bgc.restore.no2_cout = ppval(bgc.restore.no2_cs,bgc.zgrid);    
end

% N2
if ~(bgc.N2rest==0)
tdepth = -abs(Farfield.zgrid);
tconc  = Farfield.n2;
ibad = find(isnan(tconc));
tdepth(ibad) = [];
tconc(ibad) = [];
bgc.restore.n2_cs = spline(tdepth,tconc);
bgc.restore.n2_cout = ppval(bgc.restore.n2_cs,bgc.zgrid);
end

