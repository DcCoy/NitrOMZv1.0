function bgc1d_plotno2turnover(bgc)

 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 8.5 3])

 cornflowerblue = [100./255, 149./255, 237./255];
 goldenrod = [218./255,165./255,32./255];
 blue = [0, 0.4470, 0.7410];
 yellow = [0.9290, 0.6940, 0.1250];

 idx = find(bgc.no2>0.005);

 subplot(1,4,1)
 s=scatter(bgc.Data_no2(~isnan(bgc.Data_no2)), bgc.zgrid(~isnan(bgc.Data_no2)),'b');
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = [1.0 1.0 1.0];
 hold on; box on;
 plot(bgc.no2, bgc.zgrid, 'color', 'black','linewidth',3)
 ylabel('z (m)')
 xlabel('[NO_2^-] (\muM)')
 %legend('model', 'data', 'Location', 'southwest')
 
 subplot(1,4,2)
 plot(bgc.no3tono2,bgc.zgrid,'-','color',cornflowerblue,'linewidth',3)
hold on; box on;
plot(bgc.nitrox,bgc.zgrid,'-','color',goldenrod,'linewidth',3)
%ylabel('z (m)')
xlabel('rate (nM N/d)')
ylim([bgc.zbottom bgc.ztop]);
xlim([0 20]);
legend('NO_3^- red.','NO_2^- ox.', ...
     'Location', 'southeast')
 
 subplot(1,4,3)
 plot(bgc.d15no3,bgc.zgrid,'color',blue,'linewidth',3)
 hold on; box on;
 s=scatter(bgc.Data_d15no3(~isnan(bgc.Data_d15no3)), bgc.zgrid(~isnan(bgc.Data_d15no3)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = blue;
 %ylabel('z (m)')
 xlabel(insertAfter('\delta^{15}N-NO_3^- ()','(',char(8240)))
 ylim( [bgc.zbottom bgc.ztop]);

 subplot(1,4,4)
 plot(bgc.d15no2(idx),bgc.zgrid(idx),'color',yellow,'linewidth',3)
 hold on; box on;
 s=scatter(bgc.Data_d15no2(~isnan(bgc.Data_d15no2)), bgc.zgrid(~isnan(bgc.Data_d15no2)));
 s.LineWidth = 0.6;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = yellow;
 %ylabel('z (m)')
 xlabel(insertAfter('\delta^{15}N-NO_2^- ()','(',char(8240)))
 ylim( [bgc.zbottom bgc.ztop]);
 end