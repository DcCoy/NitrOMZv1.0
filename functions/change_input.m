function Tout = change_input(Tin,varargin)
% Tout = change_input(Tin,varargin)
% Change specific input values for mixed layer runs
% example :  Tout = ml_change_input(T10d_mig,'fLiMax',30);
% varargin : input pairs property/value

if mod(length(varargin),2)~=0
	error('Incorrect number of arguments - must be pairs of property/value');
end
nfields = length(varargin)/2;
varPairs = reshape(varargin,2,nfields)';

Tout = Tin;

for i=1:nfields
	thisPair = varPairs(i,:);
	thisName = thisPair{1};
	thisVal  = thisPair{2};
	% Check if property is a string
	if ~strcmp(class(thisName),'char')
		error(['Property ' num2str(i) ' should be a string']);
	end
	% Check field is in Tin
	if ~isfield(Tin,thisName)
		Tout.(thisName) = thisVal;
	else
		thisValOrig = Tin.(thisName);
		% check class of field is correct
		if class(thisVal)~=class(thisValOrig)
			error(['In property ' num2str(i) ' value is '  class(thisVal) ' kind - should be ' class(thisValOrig)]);
		end
		%substitute value
		Tout.(thisName) = thisVal;
	end
end

