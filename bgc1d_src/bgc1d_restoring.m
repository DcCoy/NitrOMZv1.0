function restoring = bgc1d_restoring(bgc,t);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Specifies restoring source (linear relaxation) for selected dissolved tracers
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% tauZvar
if bgc.tauZvar == 1
   tauh = ppval(bgc.restore.tauh_pre,bgc.zgrid);
else
   tauh = bgc.tauh;
end

% PO4
if bgc.PO4rest==0
   restoring.po4 = zeros(1,length(t.po4));
else
   restoring.po4 = (bgc.restore.po4_cout-t.po4)./tauh;
end

% NO3
if bgc.NO3rest==0
   restoring.no3 = zeros(1,length(t.no3));
else
   restoring.no3 = (bgc.restore.no3_cout-t.no3)./tauh; 
end

% O2
if bgc.O2rest==0
   restoring.o2 = zeros(1,length(t.o2));
else
   if bgc.forceanoxic == 1
      cout(find(bgc.zgrid==bgc.forceanoxic_bounds(2)):find(bgc.zgrid==bgc.forceanoxic_bounds(1)))=0;
   end
   restoring.o2 = (bgc.restore.o2_cout-t.o2)./tauh;
end

% N2O
if bgc.N2Orest==0
   restoring.n2o = zeros(1,length(t.n2o));
else
   restoring.n2o = (bgc.restore.n2o_cout-t.n2o)./tauh;
end

% NO2
if bgc.NO2rest==0
   restoring.no2 = zeros(1,length(t.no2));
else
   restoring.no2 = (bgc.restore.no2_cout-t.no2)./tauh;
end

% NH4
if bgc.NH4rest==0
   restoring.nh4 = zeros(1,length(t.nh4));
else
   restoring.nh4 = (bgc.restore.nh4_cout-t.nh4)./tauh;
end

% N2
if bgc.N2rest==0
   restoring.n2 = zeros(1,length(t.n2));
else
   cout = zeros(size(bgc.zgrid));
   restoring.n2 = (bgc.restore.n2_cout-t.n2)./tauh;
end

