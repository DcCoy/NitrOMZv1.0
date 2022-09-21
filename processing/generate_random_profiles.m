function compout = generate_random_profiles(compin,varargin) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scripts takes an existing profile compilation,
% adds some random noise
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 
%    Tracers:
%    compout = generate_random_profiles(compin,'fact',0.1);
%    Rates
%    compout = generate_random_profiles(compin,'fact',0.1, ...
%              'depth_var','depth_from_oxycline', ...
%               'var',{'no3tono2','nh4ton2o','noxton2o','nh4tono2','anammox','no2tono3'});
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default arguments
A.depth_var = 'zgrid';
A.fact = 0.1;
A.var = {'o2','no3','po4','n2o','nh4','no2','nstar'};
A = parse_pv_pairs(A,varargin);

% Process vars
varnames = A.var;
nvar = length(varnames);
tsize = size(compin.(A.depth_var));

% Initializes output structure
compout = compin;

% Seed random number generator
rng('shuffle');

% Loops through variables and adds random "noise"
for indv=1:nvar;
	tmpvar = compin.(varnames{indv});
	% Random vector
	% with values between -1 and 1
	tmpvec0 = 1-2*rand(tsize);
	% Perturb variable by relative +/- A.fact
	newvar = tmpvar .* (1 + A.fact * tmpvec0);
	% Fills in output structure
	compout.(varnames{indv}) = newvar;
	compout.([varnames{indv} '_original']) = tmpvar;
end



