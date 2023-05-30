function bgc1d_ploto2inhibition(bgc)

 figure('units','inches')
 pos = get(gcf,'pos');
 set(gcf,'pos',[pos(1) pos(2) 4 3]) % for print
 %set(gcf,'pos',[pos(1) pos(2) 10 7]) % for model evaluation

 cornflowerblue = [100./255, 149./255, 237./255];
 goldenrod = [218./255,165./255,32./255];
 teal = [71./255, 219./255, 205./255];
 gray = [128./255,128./255,128./255];

 x = linspace(0,30,300);

 plot(x, fexp(x,bgc.KO2Den1),'color',gray,'linewidth',3)
 hold on
 plot(x, fexp(x,bgc.KO2Den4),'color',cornflowerblue,'linewidth',3)
 plot(x, fexp(x,bgc.KO2Hy),'color',teal,'linewidth',3)
 plot(x, fexp(x,bgc.KO2Den2),'color',goldenrod,'linewidth',3)
 plot(x, fexp(x,bgc.KO2Den3),'color','black','linewidth',3)
 
 legend('NO_3^- \rightarrow NO_2^-', 'NO_3^- \rightarrow N_2O', ...
    'Hybrid N_2O','NO_2^- \rightarrow N_2O','N_2O \rightarrow N_2', ...
     'Location', 'northeast')
 xlabel('[O_2] (\muM)')
 ylabel('e^{-O_2/K_H^{O_2}}')

 end

