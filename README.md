# iNitrOMZ
iNitrOMZ is a nitrogen-centric biogeochemical model embeded in a below-mixed layer 1-D advection diffusion model. The model resolves a comprehensive set of processes involved in the remineralization of the sinking organic matter, starting from an imposed export flux at the base of the mixed-layer.
    
## Table of Contents

- [Updates](#updates)
- [Getting started](#getting-started)
- [Code structure](#code-structure)
- [Support](#support)
- [How to cite](#how-to-cite)

Requires MATLAB 2013 or above.

---------------------------------------
---------------------------------------
Find what you need:
[Updates](.#-Updates)
### Quick intro
### Running iNitrOMZ
### Code structure
---------------------------------------
---------------------------------------

## Updates
* 12/2019 -- Change the parameter optimization algorithm to the rather new and extremely more efficient CMAES algorithm (https://en.wikipedia.org/wiki/CMA-ES)
* 02/2020 -- Example batch jobs multi-optimization submission scripts added.
* 08/2020 -- Development branch initiated by dbianchi
* 09/2020 -- Re-organized code to store variables in bgc.Types structures, adding a "iMap" to map variables to Types, to be used for example in the "change_input" function etc.
## Getting started
#### Setting the root path
    Let's call the path to iNitrOMZ/ -- $NITROMSPATH 
      (1) Open the model initialization function
            `$NITROMSPATH/iNitrOMZ_v6.0/bgc1d_src/bgc1d_initialize.m`
            for editing and set: `bgc.root='$NITROMSPATH'`;
      (2) open the template runscript  `$NITROMSPATH/iNitrOMZv6.0/runscripts/bgc_run.m`
          for editing and set: `bgc1d_root='$NITROMSPATH/'`
#### Run the model
    Run the template script `$NITROMSPATH/iNitrOMZ_v6.0/runscripts/bgc_run.m` in MATLAB
#### Customizing the run
    Change the model defaults by modifying the initialization scripts 
    in $NITROMSPATH/bgc1d_src/ (see section on Code structure for a detailing of the 
    scripts)



## Code structure 
 #### iNitrOMZ_v6.0/runscripts/  
      Template scripts to run or optimize the model
        - bgc_run.m -- template running script
                   
 #### iNitrOMZ_v6.0/bgc1d_src/
  ##### User-customizable initialization functions
       
        - bgc1d_initialize.m -- main initialization script. The user can modify 
                                general model parameters.
        - bgc1d_initboundary.m -- the user can specify/modify boudary conditions
        - bgc1d_initbgc_params.m -- the user can specify/modify biogeochemical 
                                    parameters
        - bgc1d_initIso_params.m -- the user can specify/modify parameters related 
                                    to N isotopes
  
 ##### Core model functions 
        - bgc1d_initialize_DepParam.m -- calculates dependant model parameters
        - bgc1d_initIso_Dep_params.m -- calculates dependant model parameters 
                                        related to N isotopes
        - bgc1d_initIso_update_r15n.m -- used to update N isotopic fractions
        - bgc1d_advection_diff_opt.m -- this is the model core. This function 
                                        performs the advection and diffusion of 
                                        model tracers, applies sources and sinks,
                                        and applies restoring. Also handles model 
                                        output archiving.
        - bgc1d_sourcesink.m -- calculates sources and sinks of model tracers
        - bgc1d_restoring_initialize.m -- initializes lateral restoring forcing
        - bgc1d_restoring.m -- calculates lateral restoring of model tracers
        - not listed here: small utility functions used during intergration (e.g., n2o_yield.m)

      
 #### iNitrOMZ_v6.0/processing/ 
      Processing functions usefull for analysing the solution
        - bgc1d_postprocess.m -- processes the final archived model solution into 
                                  a user-friendly structure 
        - more not listed here ...
        
 #### iNitrOMZ_v6.0/Data/
     Forcing and validation data .
       
 #### iNitrOMZ_v6.0/restart/
      Where restart files are stored
       
 #### iNitrOMZ_v6.0/saveOut/
      Where model output is archived.
      
## Support
Contact Simon Yang or Daniele Bianchi at UCLA for support. 

## How to cite 
Please cite this repository [![DOI](https://zenodo.org/badge/236965059.svg)](https://zenodo.org/badge/latestdoi/236965059)

Manuscript reference soon to come.
