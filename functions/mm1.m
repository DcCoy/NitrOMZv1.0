function fmen1 = mmen1(var,Kvar);
% Michaelis-Menten function

% If needed turn on following line to prevent less than zero variable values
% (note that this is not necessary since the same operation is done in sms)
if (0)
   var = max(0,var);
end
fmen1 = var./(var+Kvar);

