function Y = n2o_yield_hybrid(o2,no2,bgc);

  if strcmp(bgc.n2o_yield, 'SY')    % Yang fit

  elseif strcmp(bgc.n2o_yield, 'Ji') % Ji 2018 Yield

     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   
     % % Original equation
     %    ---> N-N2o/NO2: 
     %	   Ynn2o_no2 = bgc.Ji_a./ o2 + bgc.Ji_b)/100
     %
     % %  % % % % % % % % % % % % % % % % % % % % % % % % % % % %
     % % Do some algebra:
     %   we define a1 = bgc.Ji_a/100 and b1 = bgc.Ji_b/100,
     %   ---> N-NO2/NH4:
     %	   Ynn2o_nh4 = (a1 + b1*o2)/(a1+(b1+1)*o2)
     %
     % % Now we assume that N-N2O_hx/NH4 (from Hydroxylamine)
     % % does not depend on O2, and has a value equal to
     %   ---> N-N2O_hx/NH4 = lim(o2-->inf) N-N2O/NH4:
     %      Ynn2o_nh4_hx = b1/(b1+1)
     %
     % %  % % % % % % % % % % % % % % % % % % % % % % % % % % % %
     
     % scale original params
     b1=bgc.Ji_b./100.0;
     % get yields for Hydroxylamine pathway
     Y.nn2o_hx_nh4 = b1./(b1+1);
     % get yields for hybrid pathway
     Y.nn2o_hy_nh4 = bgc.yHy .* mm1(no2,bgc.KNO2Hy) .* fexp(o2,bgc.KO2Hy); % internal NO2- source: take out .* mm1(no2,bgc.KNO2Hy)
     % calculate yield for NO2- production
     Y.no2_hx_nh4 = 1-Y.nn2o_hx_nh4- Y.nn2o_hy_nh4; % internal NO2- source: Y.no2_hx_nh4 = 1-Y.nn2o_hx_nh4- 2 .* Y.nn2o_hy_nh4;
     % Check mass balance
     if abs(Y.nn2o_hx_nh4+Y.nn2o_hy_nh4+Y.no2_hx_nh4-1) > 10^-10 
        error('Mass imbalance');
     end



  end