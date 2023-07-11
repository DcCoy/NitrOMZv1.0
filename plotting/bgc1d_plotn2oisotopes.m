function bgc1d_plotn2oisotopes(bgc, bgcsp0)

 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 8.5 3])

 cornflowerblue = [100./255, 149./255, 237./255];
 goldenrod = [218./255,165./255,32./255];
 teal = [71./255, 219./255, 205./255];
 purple = [0.4940, 0.1840, 0.5560];
 green = [0.4660, 0.6740, 0.1880];

 subplot(1,4,1)
 s=scatter(bgc.Data_n2o(~isnan(bgc.Data_n2o)), bgc.zgrid(~isnan(bgc.Data_n2o)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [1.0 1.0 1.0];
 hold on; box on;
 plot(bgc.n2o, bgc.zgrid, 'color', 'black','linewidth',3)
 ylabel('z (m)')
 xlabel('[N_2O] (\muM)')

 subplot(1,4,2)
 plot(bgc.no2ton2o,bgc.zgrid,'-','color',goldenrod,'linewidth',3)
hold on; box on;
plot(bgc.no3ton2o,bgc.zgrid,'-','color',cornflowerblue,'linewidth',3)
plot(bgc.nh4ton2o,bgc.zgrid,'-','color',teal,'linewidth',3)
%ylabel('z (m)')
xlabel('N_2O production (nM N/d)')
ylim([bgc.zbottom bgc.ztop]);
%xlim([0 20]);
%legend('NO_2^- \rightarrow N_2O','NO_3^- \rightarrow N_2O', ...
legend('NO_2^-','NO_3^-','NH_4^+', ...
     'Location', 'southeast')
 
 subplot(1,4,3)
 plot(bgcsp0.d15n2oA,bgcsp0.zgrid,'color',[.7 .7 .7],'linewidth',3)
 hold on; box on;
 plot(bgc.d15n2oA,bgc.zgrid,'color',purple,'linewidth',3)
 s=scatter(bgc.Data_d15Na(~isnan(bgc.Data_d15Na)), bgc.zgrid(~isnan(bgc.Data_d15Na)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = purple;
 %ylabel('z (m)')
 xlabel(insertAfter('\delta^{15}N-N_2O^{\alpha} ()','(',char(8240)))
 ylim( [bgc.zbottom bgc.ztop]);
 %legend(insertAfter('\delta^{15}N^{sp}=0','0',char(8240)),...
 %    insertAfter('\delta^{15}N^{sp}=22','22',char(8240)),...
 %    'Location', 'southeast')

 subplot(1,4,4)
 plot(bgc.d15n2oB,bgc.zgrid,'color',green,'linewidth',3)
 hold on; box on;
 s=scatter(bgc.Data_d15Nb(~isnan(bgc.Data_d15Nb)), bgc.zgrid(~isnan(bgc.Data_d15Nb)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = green;
% ylabel('z (m)')
 xlabel(insertAfter('\delta^{15}N-N_2O^{\beta} ()','(',char(8240)))
 ylim( [bgc.zbottom bgc.ztop]);
 end