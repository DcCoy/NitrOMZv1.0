function [norm_model norm_data] = minmax_data(constraints_model,constraints_data)
% Finds min and max of either data or model prediction, as vertical profiles
% Note: the max range (data or model) is given by: (mmax-mmin)
mmax = repmat(nanmax(constraints_data.val,[],2),1,size(constraints_model,2));
mmin = repmat(nanmin(constraints_data.val,[],2),1,size(constraints_model,2));
% Removes minimum and normalizes by range
norm_model = (constraints_model - mmin)./(mmax - mmin);
norm_data  = (constraints_data.val - mmin)./(mmax - mmin);
