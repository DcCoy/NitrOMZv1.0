function bgc1d_plot_var(bgc,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot tracers from a bgc1d_run 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default arguments
A.data = 1;
A.fact = 1;
A.var = {'o2','po4','no3','no2','nh4','n2o','nstar','n2','poc'};
A.fig = 0;
A.col = [0 0 0];
A = parse_pv_pairs(A,varargin);
A.xlabels = {'[O_2]','[PO_4^{3-}]','[NO_3^-]','[NO_2^-]','[NH_4^+]','[N_2O]','N^*','[N_2]','[POC]'};

if A.fig==0
	%hfig = figure; % for model evaluation
    hfig = figure('units','inches'); % for print
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 7 7])
else
	hfig = figure(A.fig);
end

if ~iscell(A.var)
	A.var = {A.var};
end  

% Get subplots
nvar = length(A.var);
pp = numSubplots(nvar);
sp1 = pp(1);
sp2 = pp(2);

for indv=1:nvar;
	varname = A.var{indv};
	var_plot = bgc.(varname) * A.fact; 
	if isfield(bgc,['Data_' varname])
		var_data = bgc.(['Data_' varname]);
		var_range = [min([var_plot,var_data]) max([var_plot,var_data])];
	else
		var_range = [min(var_plot) max(var_plot)];
	end
	var_span = diff(var_range);

	% Plots Variable
	subplot(sp1,sp2,indv)
	if A.data & isfield(bgc,['Data_' varname])
		s=scatter(var_data(~isnan(var_data)), bgc.zgrid(~isnan(var_data)),'b');
		s.LineWidth = 0.6;
		s.MarkerEdgeColor = 'k';
		s.MarkerFaceColor = [0.7 0.7 0.9];
		hold on; 
	end
	plot(var_plot,bgc.zgrid,'-','color',A.col,'linewidth',3)
	%title([varname])
	ylabel('z (m)')
    xlabel([A.xlabels{indv} ' (\muM)'])
	%xlabel([varname ' (\mu M)'])
	ylim([bgc.zbottom bgc.ztop]);
	xlim([var_range(1)-var_span/10 var_range(2)+var_span/10]);
	%grid on; box on;
    box on;
end
%print('../plotting/plots/bgc1d_vars','-dpng')  % temporary fix to get plots to save out
%print('-dpng',['plots/bgc1d_vars']);

