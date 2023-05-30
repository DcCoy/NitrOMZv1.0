function bgc1d_plottimesteps(bgc)

 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 5 7])

 subplot(2,2,1)
 s=scatter(bgc.Data_o2, bgc.zgrid);
 hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
 xlabel('[O_2] (\muM)')

subplot(2,2,2)
s=scatter(bgc.Data_n2o, bgc.zgrid);
hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
xlabel('[N_2O] (\muM)')

 subplot(2,2,3)
 s=scatter(bgc.Data_no3, bgc.zgrid);
 hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
xlabel('[NO_3^-] (\muM)')

 subplot(2,2,4)
 s=scatter(bgc.Data_no2, bgc.zgrid);
 hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
xlabel('[NO_2^-] (\muM)')
 
for i = 1:10
    o2 = squeeze(bgc.sol_time(i*100, 1, :));
    n2o = squeeze(bgc.sol_time(i*100, 5, :));
    no3 = squeeze(bgc.sol_time(i*100, 2, :));
    no2 = squeeze(bgc.sol_time(i*100, 7, :));
    
    color = [(1000-i*100)/1000, 0, i*100/1000];

    subplot(2,2,1)
    plot(o2, bgc.zgrid, 'color', color)
    hold on; box on;

    subplot(2,2,2)
    plot(n2o, bgc.zgrid, 'color', color)
    hold on; box on;

    subplot(2,2,3)
    plot(no3, bgc.zgrid, 'color', color)
    hold on; box on;

    subplot(2,2,4)
    plot(no2, bgc.zgrid, 'color', color)
    hold on; box on;
end


figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 5 7])

 subplot(2,2,1)
 s=scatter(bgc.Data_d15no3(~isnan(bgc.Data_d15no3)), bgc.zgrid(~isnan(bgc.Data_d15no3)));
 hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
 xlabel(insertAfter('\delta^{15}N-NO_3^- ()','(',char(8240)))

subplot(2,2,2)
s=scatter(bgc.Data_d15no2(~isnan(bgc.Data_d15no2)), bgc.zgrid(~isnan(bgc.Data_d15no2)));
hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
xlabel(insertAfter('\delta^{15}N-NO_2^- ()','(',char(8240)))

 subplot(2,2,3)
 s=scatter(bgc.Data_d15Na(~isnan(bgc.Data_d15Na)), bgc.zgrid(~isnan(bgc.Data_d15Na)));
 hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
xlabel(insertAfter('\delta^{15}N-N_2O^{\alpha} ()','(',char(8240)))

 subplot(2,2,4)
 s=scatter(bgc.Data_d15Nb(~isnan(bgc.Data_d15Nb)), bgc.zgrid(~isnan(bgc.Data_d15Nb)));
 hold on; box on;
 s.LineWidth = 0.2;
 s.MarkerEdgeColor = 'k';
 s.MarkerFaceColor = 'k';
 ylabel('z (m)')
xlabel(insertAfter('\delta^{15}N-N_2O^{\beta} ()','(',char(8240)))
 
 
for i = 1:10
    d15no3 = squeeze((bgc.sol_time(i*100, 9, :)./bgc.sol_time(i*100, 2, :)/0.003675 - 1)*1000);
    d15no2 = squeeze((bgc.sol_time(i*100, 10, :)./bgc.sol_time(i*100, 7, :)/0.003675 - 1)*1000);
    d15Na = squeeze((bgc.sol_time(i*100, 12, :)./bgc.sol_time(i*100, 5, :)/0.003675 - 1)*1000);
    d15Nb = squeeze((bgc.sol_time(i*100, 13, :)./bgc.sol_time(i*100, 5, :)/0.003675 - 1)*1000);
    
    color = [(1000-i*100)/1000, 0, i*100/1000];

    subplot(2,2,1)
    plot(d15no3, bgc.zgrid, 'color', color)
    hold on; box on;

    subplot(2,2,2)
    plot(d15no2, bgc.zgrid, 'color', color)
    hold on; box on;

    subplot(2,2,3)
    plot(d15Na, bgc.zgrid, 'color', color)
    hold on; box on;

    subplot(2,2,4)
    plot(d15Nb, bgc.zgrid, 'color', color)
    hold on; box on;
end


 end

