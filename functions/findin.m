function iivect = findin(val,vect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a value and a vector
% find the closest index in the vector to the given value
% NOTE if vector contains 2 values, it's interpreted as BOUNDARIES
% between which values should be looked for
% 
% Usage: index = findin(value,vector) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(size(vect))>2
	error(['Second input must be a 1 or 2 column array']);
end

if min(size(vect)) == 1
	for indi=1:length(val)
		dd = abs(vect-val(indi));
		[mm ii] = min(dd);
		iivect(indi) = ii(1);
	end
elseif min(size(vect)) == 2
	shortdim = find(size(vect)==2);
	longdim  = find(size(vect)~=2);
	iivect = nan(size(val));
	for indi=1:length(val)
		dd = prod(vect-val(indi),shortdim);
		ii= find(dd<=0,1,'first');
		if ~isempty(ii) 
			iivect(indi) = ii(1);
		else
			if val(indi)<min(vect(:))
				iivect(indi)=1;
			elseif val(indi)>max(vect(:))
				iivect(indi)=longdim;
			else
				disp(['WARNING: Problem encountered e.g. NaN -- Assigning NaN!']);
				iivect(indi) = nan;
			end
		end
	end
else 
	error(['Second input must be 1 or 2 column array']);
end
