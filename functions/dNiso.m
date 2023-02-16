function iso = dNiso(varargin)

A.N = nan;
A.i14N = nan;
A.i15N = nan;
A.d15N = nan;
A = parse_pv_pairs(A, varargin);

% Reference 15N/14N
Nstd1 =  0.3663/99.6337;
Nstd2 =  0.3663/100.0;

if ~isnan(A.d15N(1)) & ~isnan(A.N(1))
	iso.d15N=A.d15N; iso.N=A.N;
	tmp=(A.d15N./1000.0+1.0).*Nstd1;
	iso.i15N=(tmp.*iso.N)./(1.0+tmp);
	iso.i14N=iso.N-iso.i15N;
elseif ~isnan(A.d15N(1)) & ~isnan(A.i14N(1))
        iso.d15N=A.d15N; iso.i14N=A.i14N;
        iso.i15N=(A.d15N./1000.0+1.0).*Nstd1.*A.i14N;
	iso.N=iso.i14N+iso.i15N;
elseif ~isnan(A.N(1)) & ~isnan(A.i15N(1))
	iso.N=A.N; iso.i15N=A.i15N; iso.i14N=A.N-A.i15N;
	iso.d15N = ((A.i15N./iso.i14N)./Nstd1-1.0).*1000.0;
elseif ~isnan(A.N(1)) & ~isnan(A.i14N(1))
        iso.N=A.N; iso.i14N=A.i14N; iso.i15N=A.N-A.i14N;
        iso.d15N = ((A.i15N./iso.i14N)./Nstd1-1.0).*1000.0;
elseif ~isnan(A.i15N(1)) & ~isnan(A.i14N(1))
        iso.i14N=A.i14N; iso.N=A.i15N+A.i14N; iso.i15N=A.i15N;
        iso.d15N = ((A.i15N./iso.i14N)./Nstd1-1.0).*1000.0;	
end	

