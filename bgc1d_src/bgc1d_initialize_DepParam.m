function bgc = bgc1d_initialize_DepParam(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Calculation of parameter values dependant on those set in bgc1d_initialize
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Variables %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.varname = bgc.tracers;
bgc.nvar  = length(bgc.varname); % number of tracers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Paths %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(bgc.region, 'ETNP')
	bgc.farfield_profiles = '/Data/farfield_ETNP_gridded.mat'; % farfield concentrations
elseif strcmp(bgc.region,'ETSP')
	bgc.farfield_profiles = '/Data/farfield_ETSP_gridded.mat'; % farfield concentrations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Restart %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bgc.FromRestart == 1
	if bgc.rest_verbose
		disp(['Loading restart : ' bgc.root '/restart/' bgc.RestartFile]);
	end
	load([bgc.root,'/restart/',bgc.RestartFile]);
	bgc.rst=rst;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Vertical grid %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bgc.depth = abs(bgc.ztop-bgc.zbottom); % Water column depth (m)
bgc.nz = bgc.npt + 1; % Number of points in z
bgc.dz = (bgc.ztop - bgc.zbottom)/bgc.npt; % dz (m)
bgc.zgrid = linspace(bgc.ztop,bgc.zbottom,bgc.npt+1); % z points (m)
bgc.zgridpoc = linspace(bgc.ztop+bgc.dz/2,bgc.zbottom-bgc.dz/2,bgc.npt+2); % staggered grid for POC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Particle sinking %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If bgc.varsink = 1 then we use a martin curve remineralization profile
% We infer the depth dependant sinking speed of POC using bgc.b,  the exponent 
% of the martin curve phi = phi0*(z/z0)^b. Modify bgc.b in bgc.initialize to 
% modify the sinking profile. 
if bgc.varsink == 1					
	bgc.zpocref = bgc.zgrid(1)+bgc.dz/2	;
	bgc.a = bgc.Krem / bgc.b;
	bgc.wsinkpoc =  bgc.a .* (-bgc.zgridpoc);
	bgc.wsink =  -bgc.a .* (bgc.zgrid);    
else
	bgc.wsink = bgc.wsink_param .* ones(1,length(bgc.zgrid ));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Stochiometry %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate stoichiometry based on organic molecule: C_aH_bO_cN_dP
bgc = get_stoichiometry(bgc,bgc.stoch_a,bgc.stoch_b,bgc.stoch_c,bgc.stoch_d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Restoring %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bgc.tauZvar ~= 1
	bgc.tauh = ((bgc.Uh/bgc.Lh + 2*bgc.Kh/bgc.Lh^2)^-1)*bgc.Rh;
end 
bgc = bgc1d_restoring_initialize(bgc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Vert Diffusion %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bgc.depthvar_Kv == 1
	bgc.Kv = 0.5*(bgc.Kv_top + bgc.Kv_bot) + 0.5*(bgc.Kv_top - bgc.Kv_bot) * ...
		  tanh((bgc.zgrid - bgc.Kv_flex)/(0.5*bgc.Kv_width));
else
	bgc.Kv = bgc.Kv_param.*ones(1,length(bgc.zgrid));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Upwelling %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get upwelling profile by applying a spline to the model output w(z) profile
if bgc.depthvar_wup == 1	 
	load([bgc.root,bgc.wup_profile]);
	tdepth = -abs(vert.depth);
	if strcmp(bgc.region,'ETNP')
		tconc  = -vert.wvelZ_etnp;     
	elseif strcmp(bgc.region,'ETSP')
		tconc  = -vert.wvelZ_etsp;
	end
	ibad = find(isnan(tconc));
	tdepth(ibad) = [];
	tconc(ibad) = [];     
	wupfcn = spline(tdepth,tconc);
	bgc.wup = ppval(wupfcn,bgc.zgrid)/100.0; % convert to m/s;
else
	bgc.wup = bgc.wup_param .* ones(1,length(bgc.zgrid));
end

