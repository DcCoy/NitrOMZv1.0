function bgc = bgc1d_postprocess(bgc,Data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the ini solution from the bgc.sol structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get solution
bgc.sol = squeeze(bgc.sol_time(end,:,:));

% Get vertical derivatives
for indt=1:bgc.nvar
	bgc.(bgc.varname{indt}) = bgc.sol(indt,:);
	if bgc.flux_diag == 1
		bgc.(['adv' bgc.varname{indt}])  = bgc.sadv(indt,:);
		bgc.(['diff' bgc.varname{indt}]) = bgc.sdiff(indt,:);
		bgc.(['sms' bgc.varname{indt}])  = bgc.ssms(indt,:);
		bgc.(['rest' bgc.varname{indt}]) = bgc.srest(indt,:);
	end
	bgc.(['d' bgc.varname{indt}]) = nan(size(bgc.(bgc.varname{indt})));
	bgc.(['d' bgc.varname{indt}])(2:end-1) = (bgc.(bgc.varname{indt})(3:end) - bgc.(bgc.varname{indt})(1:end-2))/(-2*bgc.dz);
	bgc.(['d' bgc.varname{indt}])(1) = 0; 
	bgc.(['d' bgc.varname{indt}])(end) = 0;
	bgc.(['d2' bgc.varname{indt}]) = nan(size(bgc.(['d' bgc.varname{indt}])));
	bgc.(['d2' bgc.varname{indt}])(2:end-1) = (bgc.(bgc.varname{indt})(3:end) - 2 * bgc.(bgc.varname{indt})(2:end-1) + ...
		bgc.(bgc.varname{indt})(1:end-2))/(bgc.dz^2);
	bgc.(['d2' bgc.varname{indt}])(1) = 0;
	bgc.(['d2' bgc.varname{indt}])(end) = 0;
end

% Add observational data ('Data_*')
if nargin>1
	ntrData = length(bgc.varname);
	tmp = strcat('Data_', bgc.varname);
	for indt=1:ntrData
		try
			bgc.(tmp{indt}) = Data.val(indt,:);
		catch
			display(['WARNING: did not find ', tmp{indt} ' in DATA']);
		end
	end
	try
		bgc.Data_nstar = bgc.Data_no3-bgc.NCrem/bgc.PCrem*bgc.Data_po4;
	end
end

% Additional tracers:
% NSTAR - NO3 deficit versus PO4
bgc.nstar = bgc.no3 - (bgc.NCrem/bgc.PCrem) * bgc.po4;

% Get N2O advection and diffusion contributions
bgc.dwup = nan(size(bgc.wup)); bgc.dwup(1) = 0; bgc.dwup(end) = 0;
bgc.dwup(2:end-1) = (bgc.wup(3:end) - bgc.wup(1:end-2))/(-2*bgc.dz);
bgc.ssN2OAdv =  - (bgc.wup .* bgc.dn2o + bgc.dwup .* bgc.n2o);
bgc.ssN2ODiff = bgc.Kv .* bgc.d2n2o; 

% Biological sources and sinks terms
% Re-creates a "tracers" field to pass to SMS routine
for indv=1:length(bgc.tracers)
	tr.(bgc.tracers{indv}) = bgc.(bgc.tracers{indv});
end
[sms diag] = bgc1d_sourcesink(bgc,tr); 

% Converts from (uM N/s) to (nM N/d)
cnvrt = 1000*3600*24;
bgc.remox      = diag.RemOx      * cnvrt;	   % nM C/d
bgc.ammox      = diag.Ammox      * cnvrt;	   % nM n/d
bgc.anammox    = 2.0 * diag.Anammox * cnvrt;       % nM N/d : Units of N, not N2
bgc.nitrox     = diag.Nitrox     * cnvrt;	   % nM n/d
bgc.remden     = diag.RemDen     * cnvrt;	   % nM C/d
bgc.remden1    = diag.RemDen1    * cnvrt;	   % nM C/d
bgc.remden2    = diag.RemDen2    * cnvrt;	   % nM C/d
bgc.remden3    = diag.RemDen3    * cnvrt;	   % nM C/d
bgc.jn2o_ao    = diag.Jn2o_ao    * cnvrt;	   % nM N/d
bgc.jno2_ao    = diag.Jno2_ao    * cnvrt;	   % nM N/d
bgc.jn2o_prod  = 2.0 * diag.Jn2o_prod * cnvrt;     % nM N/d : Units of N, not N2O
bgc.jn2o_cons  = 2.0 * diag.Jn2o_cons * cnvrt;     % nM N/d : Units of N, not N2O
bgc.jno2_prod  = diag.Jno2_prod  * cnvrt;	   % nM n/d
bgc.jno2_cons  = diag.Jno2_cons  * cnvrt;	   % nM n/d
bgc.sms_n2o    = sms.n2o         * cnvrt;	   % nM n/d
%if bgc.RunIsotopes
% 	bgc.r15no3 = sms.r15no3;
% 	bgc.r15no2 = sms.r15no2;
% 	bgc.r15nh4 = sms.r15nh4;
% 	bgc.r15n2o = sms.r15n2o;
%	bgc.r15n2oA = bgc.i15n2oA./bgc.n2o;
%	bgc.r15n2oB = bgc.i15n2oB./bgc.n2o;
% end

% Other (for convenience)
bgc.nh4tono2 = bgc.jno2_ao; % nM N/d
bgc.no2tono3 = bgc.nitrox;  % nm N/d
% Denitrification rates
bgc.no3tono2    = bgc.NCden1 * bgc.remden1;      % nM N/d
bgc.no2ton2o    =   2 * sms.n2oind.den2 * cnvrt; % nM N/d : Units of N, not N2O
bgc.n2oton2     = - 2 * sms.n2oind.den3 * cnvrt; % nM N/d : Units of N, not NO
bgc.noxton2o    = bgc.no2ton2o;                  % nM N/d : Units of N, not N2O
% N2O formation
bgc.n2onetden   = bgc.NCden2 * bgc.remden2 - 2.0 * bgc.NCden3 * bgc.remden3; % nM N/d : Units of N, not N2O
bgc.nh4ton2o    = bgc.jn2o_ao;                                               % nM N/d : Units of N, not N2O
% Anammox fraction
bgc.AnammoxFrac = bgc.anammox ./ (bgc.anammox + bgc.NCden2.*bgc.remden2 + bgc.jn2o_prod); % non-dimensional

% Adds estimate of particle flux
% This is somewhat approximate because it's recalculated from POC and sinking speed at the tracer cells
bgc.poc_flux = -bgc.wsink .* bgc.poc*86400; % mmol C/m2/d

