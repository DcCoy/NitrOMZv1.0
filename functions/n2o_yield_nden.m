function Y = n2o_yield_nden(o2,bgc);

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
     %	---> N-NO2/NH4:	
     %	   Ynno2_nh4 = o2/(a1+(b1+1)*o2)
     %
     % % Now we assume that N-N2O_hx/NH4 (from Hydroxylamine)
     % % does not depend on O2, and has a value equal to
     %   ---> N-N2O_hx/NH4 = lim(o2-->inf) N-N2O/NH4:
     %      Ynn2o_nh4_hx = b1/(b1+1)
     %   ---> NO2_hx/NH4 = lim(o2-->inf) NO2/NH4:
     %      Ynn2o_nh4_hx = 1/(b1+1)
     %
     % % Finally, we can get N-N2O_nden/NH4 (from 
     % % nitrifier-denitrification) with
     %   ---> N-N2O_nden/NH4 = N-N2O/NH4 - N-N2O_hx/NH4
     %      Ynn2o_nh4_nden = 1/((b1+1)*(1+(b1+1)/a1)*o2))
     %   ---> N-N2O_nden/NH4 = -NO2_nden/NH4
     %      Yno2_nh4_nden = -1/((b1+1)*(1+(b1+1)/a1)*o2))
     %
     % %  % % % % % % % % % % % % % % % % % % % % % % % % % % % %
     
     % scale original params
     a1=bgc.Ji_a./100.0;b1=bgc.Ji_b./100.0;
     % get total yields
     Y.nn2o_nh4 = (a1 + b1.*o2)./(a1+(b1+1).*o2);
     Y.no2_nh4 = o2./(a1+(b1+1).*o2);
     % get yields for Hydroxylamine pathway
     Y.nn2o_hx_nh4 = b1./(b1+1);
     % get yields for nitrifier-denitrification pathway
     %Y.nn2o_nden_nh4 = 1./((b1+1).*(1+(b1+1)./a1).*o2); % CLK: I think this is a typo. Should be: 1./((b1+1).*(1+(b1+1)./a1.*o2))
     Y.nn2o_nden_nh4 =  1./((b1+1).*(1+(b1+1)./a1.*o2));
     %Y.nn2o_nden_nh4 = 1./((b1+1).*(1+(b1+1)./a1.*o2));
     Y.no2_nden_nh4 = -Y.nn2o_nden_nh4; % nh4->no2->n2o
     Y.no2_hx_nh4 = 1-Y.nn2o_hx_nh4-Y.nn2o_nden_nh4;
     % Check mass balance
     if abs(Y.nn2o_nh4+Y.no2_nh4-1) > 10^-10 
        error('Mass imbalance');
     end



  end
