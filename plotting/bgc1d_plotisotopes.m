function bgc1d_plotisotopes(bgc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bgc1d ncycle v 1.0 - Simon Yang  - October 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results from optimization by the genetic algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot modelled O2, NO3, PO4, POC, N2O vs the Data used for the optimization%%%

 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 10 7])

 idx = find(bgc.no2>0.005);
 blue = [0, 0.4470, 0.7410];
 red = [0.8500, 0.3250, 0.0980];
 yellow = [0.9290, 0.6940, 0.1250];
 purple = [0.4940, 0.1840, 0.5560];
 green = [0.4660, 0.6740, 0.1880];


 subplot(2,2,1)
 plot(bgc.d15no3,bgc.zgrid,'color',blue,'linewidth',3)
 grid on; hold on
 s=scatter(bgc.Data_d15no3(~isnan(bgc.Data_d15no3)), bgc.zgrid(~isnan(bgc.Data_d15no3)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = blue;
 ylabel('z (m)')
 xlabel(char(8240))
 title('\delta^{15}N-NO_3^-')
 ylim( [bgc.zbottom bgc.ztop]);
 legend('model', 'data', 'Location', 'southwest')

 subplot(2,2,2)
 plot(bgc.d15no2(idx),bgc.zgrid(idx),'color',red,'linewidth',3)
 grid on; hold on
 s=scatter(bgc.Data_d15no2(~isnan(bgc.Data_d15no2)), bgc.zgrid(~isnan(bgc.Data_d15no2)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = red;
 ylabel('z (m)')
 xlabel(char(8240))
 title('\delta^{15}N-NO_2^-')
 ylim( [bgc.zbottom bgc.ztop]);

 subplot(2,2,3)
 plot(bgc.d15n2oA,bgc.zgrid,'color',purple,'linewidth',3)
 grid on; hold on
 s=scatter(bgc.Data_d15Na(~isnan(bgc.Data_d15Na)), bgc.zgrid(~isnan(bgc.Data_d15Na)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = purple;
 ylabel('z (m)')
 xlabel(char(8240))
 title('\delta^{15}N-N_2O^{\alpha}')
 ylim( [bgc.zbottom bgc.ztop]);

 subplot(2,2,4)
 plot(bgc.d15n2oB,bgc.zgrid,'color',green,'linewidth',3)
 grid on; hold on
 s=scatter(bgc.Data_d15Nb(~isnan(bgc.Data_d15Nb)), bgc.zgrid(~isnan(bgc.Data_d15Nb)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = green;
 ylabel('z (m)')
 xlabel(char(8240))
 title('\delta^{15}N-N_2O^{\beta}')
 ylim( [bgc.zbottom bgc.ztop]);

 end

