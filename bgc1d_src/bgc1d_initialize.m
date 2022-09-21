function bgc = bgc1d_initialize(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Initialization of model parameters 
% Note that bgc/n-cycling parameters are defined in -- bgc1d_initbgc_params.m 
% Note that boundary condition values are defined in - bgc1d_initboundary.m
% Note that dependent parameters are updated in ------ bgc1d_initialize_DepParam.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% General %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% User specific  %%%%%%%%%
run('../bgc1d_paths_init.m');
bgc.root         = my_root; % root set in bgc1d_paths_init
bgc.RunName      = 'test_ETSP'; % set name of run
bgc.region       = 'ETSP'; % set region
bgc.wup_profile  = '/Data/vertical_CESM.mat'; % vertical velocities
bgc.Tau_profiles = '/Data/Tau_restoring.mat'; % Depth dependent Restoring timescale
bgc.visible      = 'on'; % Show figures in X window
bgc.flux_diag    = 0; % Save fluxes online
bgc.FromRestart = 0; % initialize from restart? 1 yes, 0 no
bgc.RestartFile = 'ETSP_restart.mat'; % restart file
bgc.SaveRestart = 0; %Save restart file? 1 yes, 0 no

%%%%%%%% Vertical grid %%%%%%%%%
bgc.npt = 130; % % number of mesh points for solution (for IVP)
bgc.ztop = -30; % top depth (m)
bgc.zbottom = -1330; % bottom depth (m)

%%%%% Time step / history %%%%%%
iTstep = 2;
switch iTstep
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
% 'GCN': Generalized Crank-Nicolson (currently broken April 2019)
bgc.advection = 'FTCS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Model general   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Prognostic variables %%%%%%
bgc.tracers = {'o2', 'no3','poc', 'po4', 'n2o', 'nh4', 'no2', 'n2'};
bgc.nvar_tr = length(bgc.tracers);

%%%%%%% Particle sinking %%%%%%%
bgc.varsink = 1; % if 1 then use Martin curve else, use constant sinking speed. 
if bgc.varsink == 1
	bgc.b = -0.7049; % Martin curve exponent: Pi = Phi0*(z/z0)^b
else
	bgc.wsink_param = -20/(86400); % constant speed (bgc.varsink==0)
end

%%%%%% Upwelling speed %%%%%%%%%
% Choose constant (=0) or depth-dependent (=1) upwelling velocity
% depth-dependent velocity requires a forcing file (set in bgc1d_initialize_DepParam.m)
bgc.depthvar_wup = 0; 
bgc.wup_param = 4.0 * 7.972e-8;% 1.8395e-7; % m/s  % note: 10 m/y = 3.1710e-07 m/s

%%%%%%%%%%% Diffusion %%%%%%%%%%
bgc.depthvar_Kv = 1; 
bgc.Kv_param  = 2.0 * 1.701e-5; % constant vertical diffusion coefficient in m^2/s
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
% Initialize dependent parameters. (This should be on for optimization)
bgc.depparams = 1;

% Initialize BGC/N-cycling parameters (modify in bgc1d_initbgc_params.m)
bgc = bgc1d_initbgc_params(bgc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Restoring  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up restoring timescales for far-field profiles as a crude representation
% of horizontal advective and diffusive processes.

%%%%%% On and off switches %%%%%%%%%
% Restoring switches: 1 to restore, 0 for no restoring
bgc.RestoringOff = 1;	% 1: turns restoring off for all variables
		% (supersedes all following terms, used to speedup code)
bgc.PO4rest = 0;
bgc.NO3rest = 0;
bgc.O2rest  = 0;
bgc.N2Orest = 0 ;
bgc.NH4rest = 0;
bgc.N2rest  = 0;
bgc.NO2rest = 0;

%%%%%% Z-dependent restoring timescale %%%%%%%%%
% Set to 1 for depth varying restoring timescales, 0 for constant
% bgc.tauZvar = 1 requires a forcing file set in bgc1d_initialize_DepParam.m
bgc.tauZvar = 1; 

%%%%%% Physical scalings %%%%%%
bgc.Rh = 1.0; 			% unitless scaling for sensitivity analysis. Default is 1.0
bgc.Lh = 4000.0 * 1e3;	% m - horizontal scale
% if you chose constant restoring timescales
if bgc.tauZvar == 0
	bgc.Kh = 1000;		% m2/s - horizontal diffusion
	bgc.Uh = 0.05;		% m/s - horizontal advection
end

%%%%%% Force Anoxia %%%%%%
% Option to force restoring to 0 oxygen in a certain depth range.
% Useful to force the OMZ to span a target depth range or to remove
% O2 intrusion in the OMZ while keeping restoring in the rest of the
% water column
% As usual, 1 is on and 0 is off 
bgc.forceanoxic = 0;
% Choose depth range
bgc.forceanoxic_bounds = [-350 -100]; 

% Calculate BGC/N-cycling parameters  that depend on bgc1d_initbgc_params
if bgc.depparams
	bgc = bgc1d_initialize_DepParam(bgc);
end
