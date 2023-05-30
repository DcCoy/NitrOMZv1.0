function fexp = fexp(lvar,klvar);
% Function to apply exponential O2-inhibition to anaerobic rates

% If needed turn on following line to prevent less than zero variable values
% (note that this is not necessary since the same operation is done in sms)
if (0)
	lvar = max(0,lvar);
end
fexp = exp(-lvar/klvar);

