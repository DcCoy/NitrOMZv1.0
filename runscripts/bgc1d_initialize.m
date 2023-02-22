function bgc = bgc1d_initialize(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Initialization of model parameters 
% Note that bgc/n-cycling parameters are defined in -- bgc1d_src/bgc1d_initbgc_params.m 
% Note that boundary condition values are defined in - bgc1d_src/bgc1d_initboundary.m
% Note that dependent parameters are updated in ------ bgc1d_src/bgc1d_initialize_DepParam.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% General %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root set in bgc1d_paths_init
run('../bgc1d_paths_init.m');
bgc.root        = my_root;            

%%%%%%% User specific  %%%%%%%%%
bgc.RunName      = 'spinup_ETNP'; % Set name of run
bgc.region       = 'ETNP';      % Set region ('ETSP','ETNP', or a custom region)
bgc.visible      = 'on';        % If ('on)' then show figures in window, 'off' to make invisible
bgc.flux_diag    = 0;           % If (1) then save fluxes online
bgc.FromRestart  = 0;           % If (1) then initialize from restart? (0) No
bgc.SaveRestart  = 0;           % If (1) then save restart file? (0) No
bgc.varsink      = 1;           % If (1) then use Martin curve, else use constant sinking speed. 
bgc.depthvar_wup = 0;           % Constant (0) or depth-dependent (1) upwelling velocity
bgc.depthvar_Kv  = 0; %1;       % Constant (0) or depth-dependent (1) diffusion profile %CLK: turning off depth-dependent allows us to replicate old model
bgc.iTstep       = 3;           % Constant (1) or variable (2) time-stepping
bgc.depparams    = 1;			% Initialize dependent parameters that depend on bgc1d_initbgc_params 
bgc.RestoringOff = 1;	        % If (1), turns restoring off for all variables
bgc.forceanoxic  = 0;           % If (1), force anoxia over a given depth range
bgc.tauZvar      = 1;           % If (1), use a depth-dependent restoring time-scale (requires bgc.Tau_profiles)
bgc.RunIsotopes = true; % true -> run with isotopes

%%%%%%% Data sources for wup, Tau, Restart  %%%%%%%%%
bgc.wup_profile  = '/data/vertical_CESM.mat'; % vertical velocities
bgc.Tau_profiles = '/data/Tau_restoring.mat'; % Depth dependent Restoring timescale
bgc.RestartFile  = 'ETNP_restart.mat'; % restart file

%%%%%%%% Vertical grid %%%%%%%%%
bgc.npt = 130; % % number of mesh points for solution (for IVP)
bgc.ztop = -30; % top depth (m)
bgc.zbottom = -1330; % bottom depth (m)

%%%%% Time step / history %%%%%%
switch bgc.iTstep
% Constant dt
case 1
	% Specifies # timesteps, length and hist in timesteps 
	nt = 250*365;% Simulation length in timesteps
	dt = 2.0*86400; % timestep in seconds bgc.hist =  500; 
	hist = 365*10; % save a snapshot every X timesteps
	endTimey = nt*dt/(365*86400); % end time of simulation (years)
	histTimey = hist*dt/(365*86400); % history timestep (years)
	% Creates dt and history vectors
	[dt_vec time_vec hist_time_vec hist_time_ind hist_time] = bgc1d_process_time_stepping(dt,endTimey,histTimey);
% Variable time-stepping
case 2
	dt       = [5.0 2.0 1.0 0.5 0.25 0.125]*86400/2;
	endTimey = [650 670 690 695 698 700];
	% Output time step
	histTimey = 20; % history timestep (years)
	[dt_vec time_vec hist_time_vec hist_time_ind hist_time] = bgc1d_process_time_stepping(dt,endTimey,histTimey);
case 3
	% smaller timesteps and shorter run for isotope testing
	nt = 50000;% Simulation length in timesteps
	dt = 100000; % timestep in seconds bgc.hist =  500; 
	hist = 1000; % save a snapshot every X timesteps
	endTimey = nt*dt/(365*86400); % end time of simulation (years)
	histTimey = hist*dt/(365*86400); % history timestep (years)
	% Creates dt and history vectors
	[dt_vec time_vec hist_time_vec hist_time_ind hist_time] = bgc1d_process_time_stepping(dt,endTimey,histTimey);
    bgc.dt = dt;
otherwise
	error('Timestep mode not found');
end
bgc.dt_vec = dt_vec;
bgc.hist_time_ind = hist_time_ind;
bgc.hist_time_vec = hist_time_vec;
bgc.hist_time = hist_time;
bgc.nt = length(bgc.dt_vec);
bgc.nt_hist = length(bgc.hist_time_ind);
bgc.hist_verbose = true; % prompts a message at each saving timestep
bgc.rest_verbose = true; % prompts a message when loading restart

%% Advection diffusion scheme %%
% 'FTCS': Forward in time and centered in space
bgc.advection = 'FTCS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Model general   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Prognostic variables %%%%%%
bgc.tracers = {'o2', 'no3','poc', 'po4', 'n2o', 'nh4', 'no2', 'n2'};
bgc.nvar_tr = length(bgc.tracers);
bgc.isotopes = {'i15no3', 'i15no2', 'i15nh4', 'i15n2oA', 'i15n2oB'};
bgc.nisotopes = length(bgc.isotopes);

%%%%%%% Particle sinking %%%%%%%
if bgc.varsink == 1
	bgc.b = -0.7049; % Martin curve exponent: Pi = Phi0*(z/z0)^b
else
	bgc.wsink_param = -20/(86400); % constant speed (bgc.varsink==0)
end

%%%%%% Upwelling speed %%%%%%%%%
% Depth-dependent velocity requires a forcing file (set in bgc1d_initialize_DepParam.m)
bgc.wup_param = 1.683e-7; %4.0 * 7.972e-8;% 1.8395e-7; % m/s  % note: 10 m/y = 3.1710e-07 m/s % v5.4: 1.683e-7;

%%%%%%%%%%% Diffusion %%%%%%%%%%
bgc.Kv_param  = 1.701e-5; %2.0 * 1.701e-5; % constant vertical diffusion coefficient in m^2/s % v5.4: 1.701e-5;
% For sigmoidal Kv, use the following parameters
bgc.Kv_top = 0.70 * 2.0 * 1.701e-5;
bgc.Kv_bot = 1.00 * 2.0 * 1.701e-5;
bgc.Kv_flex = -250;
bgc.Kv_width = 300;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Boundary conditions %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify in bgc1d_initboundary.m
bgc = bgc1d_initboundary(bgc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% BGC params %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize BGC/N-cycling parameters (modify in bgc1d_initbgc_params.m)
%bgc = bgc1d_initbgc_params(bgc);
bgc = bgc1d_initbgc_params_original(bgc); % replicate old version of the model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% N Isotopes params %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify in bgc1d_initIso_params.m
bgc =  bgc1d_initIso_params(bgc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Restoring  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up restoring timescales for far-field profiles as a crude representation
% of horizontal advective and diffusive processes.

%%%%%% On and off switches %%%%%%%%%
% Restoring switches: 1 to restore, 0 for no restoring. Note that RestoreOff will override these settings.
bgc.PO4rest = 0;
bgc.NO3rest = 0;
bgc.O2rest  = 0;
bgc.N2Orest = 0 ;
bgc.NH4rest = 0;
bgc.N2rest  = 0;
bgc.NO2rest = 0;
bgc.i15NO3rest = 0;
bgc.i15NO2rest = 0;
bgc.i15NH4rest = 0; 
bgc.i15N2OArest = 0;
bgc.i15N2OBrest = 0;

%%%%%% Physical scalings %%%%%%
bgc.Rh = 1.0; 			% unitless scaling for sensitivity analysis. Default is 1.0
bgc.Lh = 4000.0 * 1e3;	% m - horizontal scale
% if you chose constant restoring timescales
if bgc.tauZvar == 0
	bgc.Kh = 1000;		% m2/s - horizontal diffusion
	bgc.Uh = 0.05;		% m/s - horizontal advection
end

%%%%%% Force Anoxia %%%%%%
% Useful to force the OMZ to span a target depth range or to remove
% O2 intrusion in the OMZ while keeping restoring in the rest of the water column
% Choose depth range
bgc.forceanoxic_bounds = [-350 -100]; 

% Calculate BGC/N-cycling parameters that depend on bgc1d_initbgc_params
if bgc.depparams
	bgc = bgc1d_initialize_DepParam(bgc);
    % Calculate dependant variables relate to isotopes
	if bgc.RunIsotopes
		bgc = bgc1d_initIso_Dep_params(bgc);
	end
end
