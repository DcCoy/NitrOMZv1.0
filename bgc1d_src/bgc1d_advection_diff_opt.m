function [sol sadv sdiff ssms srest] = bgc1d_advection_diff(bgc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Advection-diffusion module
%
%    Advects, diffuses tracers and applies 
%    BGC sources minus sinks and restorings fluxes
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Initialize solutions
sol      = zeros(bgc.nt_hist,bgc.nvar,bgc.nz);
sadv     = zeros(bgc.nt_hist,bgc.nvar,bgc.nz);
sdiff    = zeros(bgc.nt_hist,bgc.nvar,bgc.nz);
ssms     = zeros(bgc.nt_hist,bgc.nvar,bgc.nz);
srest    = zeros(bgc.nt_hist,bgc.nvar,bgc.nz);
poc      = zeros(2,bgc.nz);
o2       = zeros(2,bgc.nz);
no3      = zeros(2,bgc.nz);
no2      = zeros(2,bgc.nz);
nh4      = zeros(2,bgc.nz);
n2o      = zeros(2,bgc.nz);
n2       = zeros(2,bgc.nz);
po4      = zeros(2,bgc.nz);
fpoc_out = zeros(2,bgc.nz+1);


% Run from initial conditions or from a restart file
if bgc.FromRestart == 1  
   o2(1,:)  = bgc.rst(1,:);
   no3(1,:) = bgc.rst(2,:);
   poc(1,:) = bgc.rst(3,:);
   po4(1,:) = bgc.rst(4,:);
   n2o(1,:) = bgc.rst(5,:);
   nh4(1,:) = bgc.rst(6,:);
   no2(1,:) = bgc.rst(7,:);
   n2(1,:)  = bgc.rst(8,:);  
else
   poc(1,:) = linspace(bgc.poc_flux_top/bgc.wsink(1),0.01,bgc.npt+1);
   o2(1,:)  = linspace(bgc.o2_top,bgc.o2_bot,bgc.npt+1);
   no3(1,:) = linspace(bgc.no3_top,bgc.no3_bot,bgc.npt+1);
   no2(1,:) = 10^-23;
   nh4(1,:) = 10^-23;
   n2o(1,:) = linspace(bgc.n2o_top,bgc.n2o_bot,bgc.npt+1);
   n2(1,:)  = linspace(bgc.n2_top,bgc.n2_bot,bgc.npt+1);
   po4(1,:) = linspace(bgc.po4_top,bgc.po4_bot,bgc.npt+1);
end

% Dump tracers in a structure "tr" - one by one (avoids eval)
tr.o2  = o2(1,:);
tr.no3 = no3(1,:);
tr.poc = poc(1,:);
tr.po4 = po4(1,:);
tr.n2o = n2o(1,:);
tr.nh4 = nh4(1,:);
tr.no2 = no2(1,:);
tr.n2  = n2(1,:);

% Get initial SMS
sms =  bgc1d_sourcesink(bgc,tr);

% Initialize particulate flux at the top
fpoc_out(1,1) = bgc.poc_flux_top;

%  Update steady state POC sinking flux
for indz=1:bgc.nz
   % Updates POC, by mass conservation: remin = divergence of flux
   fpoc_out(1,indz+1) = fpoc_out(1,indz)/(1.0 + bgc.dz*sms.kpoc(indz)/bgc.wsink(indz));
   poc(1,indz)        = (fpoc_out(1,indz)-fpoc_out(1,indz+1))/(bgc.dz*sms.kpoc(indz));
end

% Start time-stepping
for indt=1:bgc.nt
   % Gets current timestep
   dt = bgc.dt_vec(indt);

   % For advection velocity and diffusion coefficient fixed in time, calculate here
   % terms for the numerical advection-diffusion solver. For time-dependent w and Kv
   % move these terms inside the time loop
   alpha = bgc.wup(2:end-1) * dt / (2*bgc.dz);
   beta  = - dt / (2*bgc.dz) * (bgc.wup(1:end-2) - bgc.wup(3:end));
   gamma = bgc.Kv(2:end-1) * dt / (bgc.dz)^2;
   delta =   dt / (4*bgc.dz) * (bgc.Kv(1:end-2) - bgc.Kv(3:end));

   % Integration coefficients for the tracer at k,k+1,k-1 vertical levels:
   coeff1 = 1 + beta - 2*gamma;
   coeff2 =     alpha +  gamma - delta;
   coeff3 =   - alpha +  gamma + delta;

   % Now calculate Explicit tracer concentrations
   % Top boundary conditions
   o2(2,1)  = bgc.o2_top;
   no3(2,1) = bgc.no3_top;
   no2(2,1) = bgc.no2_top;
   nh4(2,1) = bgc.nh4_top;
   n2o(2,1) = bgc.n2o_top;
   n2(2,1)  = bgc.n2_top;
   po4(2,1) = bgc.po4_top;
   % Bottom boundary conditions
   o2(2,end)  = bgc.o2_bot;
   no3(2,end) = bgc.no3_bot;
   no2(2,end) = bgc.no2_bot;
   nh4(2,end) = bgc.nh4_bot;
   n2o(2,end) = bgc.n2o_bot;
   n2(2,end)  = bgc.n2_bot;
   po4(2,end) = bgc.po4_bot;

   % advection and diffusion
   % The code below more compactly and efficiently solve the following equation for all tracers:
   % o2(2,2:end-1) = o2(1,2:end-1) -bgc.wup(2:end-1).*bgc.dt./(2.*-bgc.dz) .* (o2(1,3:end)-o2(1,1:end-2)) - ...
   %                 o2(1,2:end-1) .* bgc.dt./(2.*-bgc.dz) .* (bgc.wup(3:end)-bgc.wup(1:end-2)) +  ...
   %                 bgc.Kv(2:end-1) .* bgc.dt./(bgc.dz)^2 .* (o2(1,3:end) - 2 .* o2(1,2:end-1) + o2(1,1:end-2)); 

   % Explicitly goes through tracers
   o2(2,2:end-1)  = o2(1,2:end-1)  .* coeff1 + o2(1,3:end)  .* coeff2 + o2(1,1:end-2)  .* coeff3;
   no3(2,2:end-1) = no3(1,2:end-1) .* coeff1 + no3(1,3:end) .* coeff2 + no3(1,1:end-2) .* coeff3;
   po4(2,2:end-1) = po4(1,2:end-1) .* coeff1 + po4(1,3:end) .* coeff2 + po4(1,1:end-2) .* coeff3;
   n2o(2,2:end-1) = n2o(1,2:end-1) .* coeff1 + n2o(1,3:end) .* coeff2 + n2o(1,1:end-2) .* coeff3;
   nh4(2,2:end-1) = nh4(1,2:end-1) .* coeff1 + nh4(1,3:end) .* coeff2 + nh4(1,1:end-2) .* coeff3;
   no2(2,2:end-1) = no2(1,2:end-1) .* coeff1 + no2(1,3:end) .* coeff2 + no2(1,1:end-2) .* coeff3;
   n2(2,2:end-1)  = n2(1,2:end-1)  .* coeff1 + n2(1,3:end)  .* coeff2 + n2(1,1:end-2)  .* coeff3;

   % Get sources minus Sinks   
   % dump tracers in a structure "tr" - one by one (avoids eval)
   tr.o2  = o2(1,:);
   tr.no3 = no3(1,:);
   tr.poc = poc(1,:);
   tr.po4 = po4(1,:);
   tr.n2o = n2o(1,:);
   tr.nh4 = nh4(1,:);
   tr.no2 = no2(1,:);
   tr.n2  = n2(1,:);

   % Calculate SMS
   sms =  bgc1d_sourcesink(bgc,tr);

   % Update steady state POC sinking flux
   fpoc_out(2,1) = bgc.poc_flux_top;

   for indz=1:bgc.nz
      fpoc_out(2,indz+1) = fpoc_out(2,indz)/(1.0 + bgc.dz*sms.kpoc(indz)/bgc.wsink(indz));
   end       
   % Use array calculation to perform POC flux update (faster)
   poc(2,:) = (fpoc_out(2,1:bgc.nz)-fpoc_out(2,2:bgc.nz+1))./(bgc.dz*sms.kpoc);

   % Do sources minus sinks
   o2(2,2:end-1)  = o2(2,2:end-1)  + sms.o2(2:end-1)  * dt;
   no3(2,2:end-1) = no3(2,2:end-1) + sms.no3(2:end-1) * dt;
   no2(2,2:end-1) = no2(2,2:end-1) + sms.no2(2:end-1) * dt;
   nh4(2,2:end-1) = nh4(2,2:end-1) + sms.nh4(2:end-1) * dt;
   n2o(2,2:end-1) = n2o(2,2:end-1) + sms.n2o(2:end-1) * dt;
   n2(2,2:end-1)  = n2(2,2:end-1)  + sms.n2(2:end-1)  * dt;
   po4(2,2:end-1) = po4(2,2:end-1) + sms.po4(2:end-1) * dt;

   % Do restoring    
   if bgc.RestoringOff~=1
      restoring = bgc1d_restoring(bgc,tr);
      o2(2,2:end-1)  = o2(2,2:end-1)  + restoring.o2(2:end-1)  * dt;
      no3(2,2:end-1) = no3(2,2:end-1) + restoring.no3(2:end-1) * dt;
      no2(2,2:end-1) = no2(2,2:end-1) + restoring.no2(2:end-1) * dt;
      nh4(2,2:end-1) = nh4(2,2:end-1) + restoring.nh4(2:end-1) * dt;
      n2o(2,2:end-1) = n2o(2,2:end-1) + restoring.n2o(2:end-1) * dt;
      n2(2,2:end-1)  = n2(2,2:end-1)  + restoring.n2(2:end-1)  * dt;
      po4(2,2:end-1) = po4(2,2:end-1) + restoring.po4(2:end-1) * dt;    
   end

   % old equals new  
   o2(1,:)  = o2(2,:);
   no3(1,:) = no3(2,:);
   no2(1,:) = no2(2,:);
   nh4(1,:) = nh4(2,:);
   n2o(1,:) = n2o(2,:);
   n2(1,:)  = n2(2,:);
   po4(1,:) = po4(2,:);
   poc(1,:) = poc(2,:);  

   % Save history files and diagnostics  
   if any(bgc.hist_time_ind == indt)
      iout = find(bgc.hist_time_ind == indt);
      if bgc.hist_verbose
         disp(['Saving step #' num2str(iout) '/' num2str(bgc.nt_hist)]);
      end
      % Saving tracer field
      sol(iout,1,:) = o2(1,:);
      sol(iout,2,:) = no3(1,:);
      sol(iout,3,:) = poc(1,:);
      sol(iout,4,:) = po4(1,:);
      sol(iout,5,:) = n2o(1,:);
      sol(iout,6,:) = nh4(1,:);
      sol(iout,7,:) = no2(1,:);
      sol(iout,8,:) = n2(1,:);

      % Save fluxes (bgc.flux_diag == 1)
      if bgc.flux_diag == 1

         % Save advection terms
         sadv(iout,1,2:end-1) = alpha .* (o2(1,3:end)  - o2(1,1:end-2))  + beta .* o2(1,2:end-1);
         sadv(iout,2,2:end-1) = alpha .* (no3(1,3:end) - no3(1,1:end-2)) + beta .* no3(1,2:end-1);
         sadv(iout,3,2:end-1) = alpha .* (po4(1,3:end) - po4(1,1:end-2)) + beta .* po4(1,2:end-1);
         sadv(iout,4,2:end-1) = alpha .* (n2o(1,3:end) - n2o(1,1:end-2)) + beta .* n2o(1,2:end-1);
         sadv(iout,6,2:end-1) = alpha .* (nh4(1,3:end) - nh4(1,1:end-2)) + beta .* nh4(1,2:end-1);
         sadv(iout,7,2:end-1) = alpha .* (no2(1,3:end) - no2(1,1:end-2)) + beta .* no2(1,2:end-1);
         sadv(iout,8,2:end-1) = alpha .* (n2(1,3:end)  - n2(1,1:end-2))  + beta .* n2(1,2:end-1);

         % Save diffusion terms
         sdiff(iout,1,2:end-1) = gamma .* (o2(1,3:end)  - 2 * o2(1,2:end-1)  + o2(1,1:end-2));
         sdiff(iout,2,2:end-1) = gamma .* (no3(1,3:end) - 2 * no3(1,2:end-1) + no3(1,1:end-2));
         sdiff(iout,3,2:end-1) = gamma .* (po4(1,3:end) - 2 * po4(1,2:end-1) + po4(1,1:end-2));
         sdiff(iout,4,2:end-1) = gamma .* (n2o(1,3:end) - 2 * n2o(1,2:end-1) + n2o(1,1:end-2));
         sdiff(iout,6,2:end-1) = gamma .* (nh4(1,3:end) - 2 * nh4(1,2:end-1) + nh4(1,1:end-2));
         sdiff(iout,7,2:end-1) = gamma .* (no2(1,3:end) - 2 * no2(1,2:end-1) + no2(1,1:end-2));
         sdiff(iout,8,2:end-1) = gamma .* (n2(1,3:end)  - 2 * n2(1,2:end-1)  + n2(1,1:end-2));

         % Save SMS term
         ssms(iout,1,2:end-1) = sms.o2(2:end-1);
         ssms(iout,2,2:end-1) = sms.no3(2:end-1);
         ssms(iout,4,2:end-1) = sms.po4(2:end-1);
         ssms(iout,5,2:end-1) = sms.n2o(2:end-1);
         ssms(iout,6,2:end-1) = sms.nh4(2:end-1);
         ssms(iout,7,2:end-1) = sms.no2(2:end-1);
         ssms(iout,8,2:end-1) = sms.n2(2:end-1);

         % Save restoring term
         srest(iout,1,2:end-1) = restoring.o2(2:end-1);
         srest(iout,2,2:end-1) = restoring.no3(2:end-1);
         srest(iout,4,2:end-1) = restoring.po4(2:end-1);
         srest(iout,5,2:end-1) = restoring.n2o(2:end-1);
         srest(iout,6,2:end-1) = restoring.nh4(2:end-1);
         srest(iout,7,2:end-1) = restoring.no2(2:end-1);
         srest(iout,8,2:end-1) = restoring.n2(2:end-1);
      end
   end 
end  

% Save restart (bgc.SaveRestart == 1)
if bgc.SaveRestart == 1
   rst = squeeze(sol(end,:,:));
   endtime = num2str(bgc.hist_time_vec(end)/3600/24/365,'%5.1f');
   save([bgc.root, '/restart/', bgc.RunName,'_restart_',endtime,'.mat'],'rst');
end   

