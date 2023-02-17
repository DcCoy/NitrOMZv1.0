function bgc = bgc1d_run(varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template iNitrOMZ runscript 
% Usage:
%	- bgc = bgc1d_run(varargin)
%
% Inputs:
% - iPlot = (1) to plot input 
% - ParNames = Pass parameter names for update
% - ParVal   = Corresponding parameter values for 'ParNames'
%
%
% Customize your model run in bgc.root/UserParams/
%   % General model set-up		 -- bgc1d_initialize.m
%   % Boundary conditions        -- bgc1d_initboundary.m
%   % BGC/N-cycling params       -- bgc1d_initbgc_params.m
%   % N-isotopes-cycling params  -- bgc1d_initIso_params.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add root path % CLK: moved this up so it's first
run('../bgc1d_paths_init.m');
addpath(genpath(my_root));

% Process inputs (varargin)
A.iPlot = 0;			% To plot output
A.ParNames = {};		% Pass parameter names that need to be modified from default values
A.ParVal = [];			% Pass parameter values, corresponding to ParNames
A = parse_pv_pairs(A, varargin);

% Initialize the model
clear bgc;
bgc = struct;
bgc = bgc1d_initialize(bgc); 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% In case parameters are specified as inputs, e.g. by passing the
% results of an Optimization then substitutes parameters
% ParNames = {'par1';'par2';...}
% ParVal = [val1; val2; ...];
% % % % % % % % % % % % % % % % % % % % 
% Substitute any parameters from their values in bgc1d_src/bgc1d_initbgc_params.m
if ~isempty(A.ParNames) 
	for indp=1:length(A.ParNames)
	   bgc = change_input(bgc,A.ParNames{indp},A.ParVal(indp));
	end
	% Updates BGC/N-cycling parameters  that depend on bgc1d_initbgc_params
	if bgc.depparams
	   bgc = bgc1d_initialize_DepParam(bgc);
	end
end

% Run the model 
% % % % % % % % % % % % % % % % % % % % % % % % 
%     bgc.sol_time is a TxVxZ matrix where T is archive times
%     V is the number of tracers and Z in the number of 
%	  model vertical levels
%     note that the model saves in order:
%	  (1) o2 (2) no3 (3) poc (4) po4 (5) n2o (6) nh4 (7) no2 (8) n2
% % % % % % % % % % % % % % % % % % % % % % % % 
tic;
[bgc.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc);
bgc.RunTime = toc;
disp(['Runtime : ' num2str(bgc.RunTime)]);

% Process observations to validate the model solution
Tracer.name = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'};
if strcmp(bgc.region,'ETNP')
	load([bgc.root,'/data/compilation_ETNP_gridded.mat']);
	Data = proc_data(bgc,compilation_ETNP_gridded,Tracer.name);
elseif strcmp(bgc.region,'ETSP')
	load([bgc.root,'/data/compilation_ETSP_gridded.mat']);
	Data = proc_data(bgc,compilation_ETSP_gridded,Tracer.name);
end

% Process model output for analysis (gathers tracers and diagnostics into the bgc structure)
bgc = bgc1d_postprocess(bgc, Data);
if bgc.RunIsotopes
    bgc = bgc1d_pprocess_isotopes(bgc);
end
if (A.iPlot)
	bgc1d_plot_vars(bgc); % CLK: change from "bgc1d_plot" to "bgc1d_plot_vars"
    bgc1d_plot_rates(bgc); % CLK: add plot rates
end

