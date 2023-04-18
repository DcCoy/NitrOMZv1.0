function bgc = get_stochiometry(bgc,a,b,c,d)
% Calculates the stochiometry of various redox reactions given the
% composition of organic matter i.e., (C_a)(H_b)(O_c)(N_d)(P)

% Number of electrons required to oxidize (C_a)(H_b)(O_c)(N_d)(P) to CO2 +
% NH3 + H2O + H3PO4
Corg_e = 4*a+b-2*c-3*d+5;
% Number of electrons required to reduce O2 to H2O
O2toH2O_e = 4;
% number of electrons required to reduce HNO3 to HNO2
HNO3toHNO2_e = 2;
% Number of electrons required to reduce HNO2 to N2O ---> Need to check b/c # of N atoms changes
HNO2toN2O_e = 4; % CLK: 4 electrons because it's 2 e- per HNO2, times 2 HNO2
%numbers of electrons required to reduce HNO3 to N2O
HNO3toN2O_e = 8;
% Number of electrons required to reduce N2O to N2
N2OtoN2_e = 2;

% Get conversion factors
bgc.OCrem = (Corg_e / O2toH2O_e) / a;          % Oxygen to Carbon remineralization ratio (molO2/molC)
bgc.NCrem = d / a;                             % Ammonium to Carbon remineralization ratio (molNH4/molC)
bgc.PCrem = 1 / a;                             % Phosphate to Carbon denitrification ratio (molNO3/molC)
bgc.NCden1 = (Corg_e / HNO3toHNO2_e) / a;      % Nitrate to Carbon ratio during nitrate reduction (molNO3/molC)
bgc.NCden2 = (Corg_e / HNO2toN2O_e) / a;       % Nitrite to Carbon ratio during nitrite reduction to N2O (molNO2/molC)
bgc.NCden3 = (Corg_e / N2OtoN2_e) / a;         % N2O to Carbon ratio during N2O reduction (molN2O/molC)
bgc.NCden4 = (Corg_e / HNO3toN2O_e) / a;         % Nitrate to Carbon ratio during nitrate reduction to N2O (molNO3/molC)
bgc.PCden1 = 1 / a;
bgc.PCden2 = 1 / a;
bgc.PCden3 = 1 / a;

end