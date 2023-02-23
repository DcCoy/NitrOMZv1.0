function bgc_suite = bgc1d_run_suite(ParName,ParVals,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template iNitrOMZ suite runscript 
% NOTE: Requires Parallel Computing Toolbox
% Usage:
%	- bgc_suite = bgc1d_run_suite(ParName,ParVals)
%
% Inputs:
%	- ParName  = Parameter to vary
%	- ParVals  = Array of parameter values

% Optional inputs (varargin):
%	- parallel = Number of parallel pool workers 
%
% Outputs:
%	- bgc_suite = cell array of size(ParVals) with bgc1d_run outputs
%
% Example:
%	- bgc_suite = bgc1d_run_suite('KDen1',linspace(0.01,1,10));
%	- bgc_suite = bgc1d_run_suite('KDen1',linspace(0.01,1,10),'parallel',10);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add root path
run('../bgc1d_paths_init.m');
addpath(genpath(my_root));

% Process inputs (varargin)
A.parallel = [];          
A = parse_pv_pairs(A, varargin);

% Initialize the model
clear bgc;
bgc = struct;
bgc = bgc1d_initialize(bgc); 

% Initialize suite output, and change parameter values 
for i = 1:length(ParVals)
	bgc_suite{i} = bgc;
	bgc_suite{i} = change_input(bgc_suite{i},ParName,ParVals(i));
	% Updates BGC/N-cycling parameters  that depend on bgc1d_initbgc_params
	if bgc_suite{i}.depparams
	   bgc_suite{i} = bgc1d_initialize_DepParam(bgc_suite{i});
	end
end

% Check for parallel switch
if ~isempty(A.parallel)
	% Create parallel pool and run
	delete(gcp('nocreate'));
	parpool(A.parallel);
	parfor i = 1:length(ParVals);
		tic;
		[bgc_suite{i}.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc_suite{i});
		bgc_suite{i}.RunTime = toc;	
	end
	delete(gcp);
else
	% Run individually
	for i = 1:length(ParVals)
		tic;
        [bgc_suite{i}.sol_time, ~, ~, ~, ~] = bgc1d_advection_diff_opt(bgc_suite{i});
        bgc_suite{i}.RunTime = toc;
    end
end

% Process observations to validate the model solution
Tracer.name = {'o2' 'no3' 'poc' 'po4' 'n2o' 'nh4' 'no2' 'n2'};
if strcmp(bgc_suite{1}.region,'ETNP')
    load([bgc_suite{1}.root,'/data/compilation_ETNP_gridded.mat']);
	for i = 1:length(ParVals)
		Data{i} = proc_data(bgc_suite{i},compilation_ETNP_gridded,Tracer.name);
		bgc_suite{i} = bgc1d_postprocess(bgc_suite{i}, Data{i});
	end
elseif strcmp(bgc.region,'ETSP')
    load([bgc_suite{1}.root,'/data/compilation_ETSP_gridded.mat']);
	for i = 1:length(ParVals)
		Data{i} = proc_data(bgc,compilation_ETSP_gridded,Tracer.name);
		bgc_suite{i} = bgc1d_postprocess(bgc_suite{i}, Data{i});
	end
end	

