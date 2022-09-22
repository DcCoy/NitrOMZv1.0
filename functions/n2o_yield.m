function Y = n2o_yield(o2,bgc);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   
% Get N2O/NO2 yields from ammonium oxidation
% Based on parameterization from Nevison et al., 2003
%
% Y_n2o/Y_no2 = (a/o2 + b) * 0.01
% and
% Y_n2o + Y_no2 = 1
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   

switch bgc.n2o_yield
case {'Ji','Yang'}
	% Scale original params
	a1 = bgc.Ji_a ./ 100.0;
	b1 = bgc.Ji_b ./ 100.0;
	% Get total yields
	Y.n2o = (a1+b1.*o2)./(a1+(b1+1).*o2);
	Y.no2 = o2./(a1+(b1+1).* o2);
otherwise
	error(['N2O yield case : ' bgc.n2o_yield ' not found']);
end
