# iNitrOMZ
iNitrOMZ is a nitrogen-centric biogeochemical model embeded in a below-mixed layer 1-D advection diffusion model. The model resolves a comprehensive set of processes involved in the remineralization of the sinking organic matter, starting from an imposed export flux at the base of the mixed-layer.
    
## Table of Contents

- [Updates](#updates)
- [Getting started](#getting-started)
- [Code structure](#code-structure)
- [Support](#support)
- [How to cite](#how-to-cite)

Requires MATLAB 2013 or above.

## Updates
* 09/2022 -- First commit of NitrOMZv1.0 

## Getting started
#### Update settings in runscripts/bgc1d_initialize.m
	Set run-specific settings, such as:
		RunName      = Set the name of the run
		region       = Set region (i.e. 'ETSP')
		FromRestart  = Switch to initialize run from a restart file
		SaveRestart  = Switch to save output as a restart file
		varsink      = Switch to use a Martin Curve for POC sinking speed, or a constant speed
		depthvar_wup = Switch to use a constant or depth-dependent vertical upwelling velocity
		depthvar_Kv  = Switch to use a constant or depth-dependent vertical diffusion profile
		iTstep       = Switch to use a constant or variable time-step
		depparams    = Switch to initialize depth-dependent parameters
		RestoringOff = Switch to remove all restoring terms
		forceanoxic  = Switch to force anoxia over a specific depth range
	The user can also toggles settings for:
		The vertical grid: (top, bottom, and dz)
		Time-stepping: (dt, and rate of history output)
		POC sinking speed: (Martin 'b' value, or constant speed value)
		Diffusion: (constant value, or the depth-dependent parameters that control the Kv shape)
		Restoring: (which tracers will have restoring turned on, and restoring time-scale parameters)
		Anoxia: (depth bounds of forced anoxia)
#### Update boundary conditions in runscripts/bgc1d_initboundary.m
	Modify any top/bottom boundary conditions for bgc.region 
#### Update stoichiometric ratios and bgc parameters in runscripts/bgc1d_initbgc_params.m
	Modify the chemical form of organic matter (currently Anderson & Sarmiento 1994)
	Update maximum rates, half-saturation constants for oxidants/reductants, and O2 inhibition parameters
	Update parameters for calculation of N2O/NO2 yields from ammonium oxidation 
#### Run the model using runscripts/bgc1d_run
    Run the template script in MATLAB
		bgc = bgc1d_run;

## Code structure 
#### runscripts/  
	Folder where NitrOMZ runscripts and run settings files are stored
		bgc1d_run.m            -- script to run the model
		bgc1d_initialize.m     -- script to toggle run settings
		bgc1d_initboundary.m   -- script to toggle region-specific top/bottom boundary conditions
		bgc1d_initbgc_params.m -- script to toggle nitrogen cycle parameters
                   
#### bgc1d_run/  
	Folder where core model functions are stored
		bgc1d_advection_diff_opt.m    -- NitrOMZ advection/diffusion module
		bgc1d_initialize_DepParam.m   -- Calculates dependent model parameter values
		bgc1d_process_time_stepping.m -- Function to calculate time-stepping variables
		bgc1d_restoring.m             -- Applies lateral restoring forcing for selected (or all) tracers
		bgc1d_restoring_initialize.m  -- Initializes lateral restoring forcing
		bgc1d_sourcesink.m            -- Calculates sources and sinks for model tracers

#### processing/ 
	Folder where post-processing functions are stored
		bgc1d_postprocess.m -- Process the final archived model soltution into a user-friendly structure
        
#### data/
	Region-specific forcing and validation data
		compilation_ETNP_gridded.mat -- ETNP observed tracers, for validation 
		compilation_ETSP_gridded.mat -- ETSP observed tracers, for validation
		farfield_ETNP_gridded.mat    -- ETNP lateral restoring tracer concentrations 
		farfield_ETNP_gridded.mat    -- ETSP lateral restoring tracer concentrations
		comprates_ETSP.mat           -- ETSP observed N transformation rates, for validation
		Tau_restoring.mat            -- Region-specific restoring time-scales
		vertical_CESM.mat            -- Region-specific vertical upwelling data

#### restart/
	Folder to store restart files

#### plotting/
	Folder containing post-run plotting scripts
		bgc1d_plot_rates -- Plots the final N-cycle rate profiles
		bgc1d_plot_vars  -- Plots the final N-cycle tracer profiles
      
#### functions/
	Folder containing additional functions needed to run the model

## Support
Contact Daniel McCoy or Daniele Bianchi at UCLA for support. 

## How to cite 
Please cite this repository [![DOI](https://zenodo.org/badge/236965059.svg)](https://zenodo.org/badge/latestdoi/236965059)

Manuscript reference soon to come.
