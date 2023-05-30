function bgc1d_plot_rate(bgc,rates,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot rates from a bgc1d_run 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default arguments
A.data = 1;
A.fact = 1;
A.oxy_threshold = 1.0;	% o2 threshold ofr oxycline
A.mode = 'oxycline'; % 'oxycline' references depth to oxycline depth (default)
%A.var = {'ammox','nitrox','nh4ton2o','n2onetden','no3tono2','no2ton2o','n2oton2','anammox'};
A.var = {'ammox','nitrox','nh4ton2o','no3tono2','no3ton2o','no2ton2o','n2oton2','anammox'};
A.xlabels = {'NH_4^+ \rightarrow NO_2^-','NO_2^- \rightarrow NO_3^-',...
    'NH_4^+ \rightarrow N_2O', 'NO_3^- \rightarrow NO_2^-',...
    'NO_3^- \rightarrow N_2O','NO_2^- \rightarrow N_2O',...
    'N_2O \rightarrow N_2','anammox'};
A.fig = 0;
A.col = [0 0 0];
A = parse_pv_pairs(A,varargin);

if A.fig==0
	%hfig = figure; % for model evaluation
    hfig = figure('units','inches'); % for print
    pos = get(gcf,'pos');
    set(gcf,'pos',[pos(1) pos(2) 8.5 7])
else
	hfig = figure(A.fig);
end

if ~iscell(A.var)
	A.var = {A.var};
end  
A.var = lower(A.var);

% Process observational rates -- here assumes the coordinate is "depth_from_oxycline"
nan_vect = nan(1,bgc.nz);
% Adds in the rates from the data
if strcmp(A.mode,'oxycline')
	% by referring them to the oxycline them  
	ind_o2 = find(strcmp(bgc.varname,'o2'));
	depthox = bgc1d_detect_oxycline(bgc.sol(ind_o2,:),bgc,A.oxy_threshold);
	if ~isnan(depthox(1))
		for indv=1:length(A.var)
			try
				grid_data = bgc1d_griddata(bgc.(A.var{indv}),bgc.depth_from_oxycline+depthox(1),bgc);
				bgc.(['Data_' A.var{indv}]) = grid_data;
			catch
				bgc.(['Data_' A.var{indv}]) = nan_vect;
			end
		end
	end
else
	for indv=1:length(A.var)
		bgc.(['Data_' A.var{indv}]) = nan_vect;
	end
end

% Get subplots
nvar = length(A.var);
pp = numSubplots(nvar);
sp1 = pp(1);
sp2 = pp(2);

for indv=1:nvar;
	varname = A.var{indv};
	var_plot = bgc.(varname) .* A.fact; 
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
    xlabel([A.xlabels{indv} ' (nM N/d)'])
	%xlabel([varname ' nM N/d'])
	ylim([bgc.zbottom bgc.ztop]);
	xlim([var_range(1)-var_span/10 var_range(2)+var_span/10]);
	%grid on; box on;
    box on;
end
%print('../plotting/plots/bgc1d_rates','-dpng'); % temporary fix to get plots to save out
%print('-dpng',['plots/bgc1d_rates']);


