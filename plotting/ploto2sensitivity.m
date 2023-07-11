o2vals = linspace(33,200,10);

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
 
 for i = 1:length(o2vals)
    o2 = o2vals(i);
    savestr = insertAfter("../saveout/o2sensitivity.mat","sensitivity",string(o2));
    color = [(o2-33)/200, 0, (233-o2)/233];
    disp(savestr)
    bgc = load(savestr).bgc;

    subplot(2,2,1)
    plot(bgc.o2, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,2)
    plot(bgc.n2o, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,3)
    plot(bgc.no3, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,4)
    plot(bgc.no2, bgc.zgrid, 'color', color,'linewidth',3)
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
 
 
 for i = 1:length(o2vals)
    o2 = o2vals(i);
    savestr = insertAfter("../saveout/o2sensitivity.mat","sensitivity",string(o2));
    color = [(o2-33)/200, 0, (233-o2)/233];
    disp(savestr)
    bgc = load(savestr).bgc;

    subplot(2,2,1)
    plot(bgc.d15no3, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,2)
    plot(bgc.d15no2, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,3)
    plot(bgc.d15n2oA, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,4)
    plot(bgc.d15n2oB, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;
 end

  figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 5 7])

 subplot(2,2,1)
 hold on; box on;
 ylabel('z (m)')
 xlabel('NO_2^- \rightarrow N_2O ')

subplot(2,2,2)
hold on; box on;
 ylabel('z (m)')
xlabel('NO_3^- \rightarrow N_2O')

 subplot(2,2,3)
 hold on; box on;
 ylabel('z (m)')
xlabel('NH_4^+ \rightarrow N_2O')

 subplot(2,2,4)
 hold on; box on;
 ylabel('z (m)')
xlabel('N_2O \rightarrow N_2')
 
 for i = 1:length(o2vals)
    o2 = o2vals(i);
    savestr = insertAfter("../saveout/o2sensitivity.mat","sensitivity",string(o2));
    color = [(o2-33)/200, 0, (233-o2)/233];
    disp(savestr)
    bgc = load(savestr).bgc;

    subplot(2,2,1)
    plot(bgc.no2ton2o, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,2)
    plot(bgc.no3ton2o, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,3)
    plot(bgc.nh4ton2o, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;

    subplot(2,2,4)
    plot(bgc.n2oton2, bgc.zgrid, 'color', color,'linewidth',3)
    hold on; box on;
end
