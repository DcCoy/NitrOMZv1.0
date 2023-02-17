function bgc = bgc1d_pprocess_isotopes(bgc)

for indt=1:bgc.nisotopes
    bgc.(bgc.isotopes{indt}) = bgc.sol(indt+8,:);
end

bgc.d15no3 = (bgc.i15no3(:)./bgc.no3(:)/0.003675 - 1)*1000;
bgc.d15no2 = (bgc.i15no2(:)./bgc.no2(:)/0.003675 - 1)*1000;
bgc.d15nh4 = (bgc.i15nh4(:)./bgc.nh4(:)/0.003675 - 1)*1000;

idx = find(bgc.no3 < bgc.IsoThreshold | bgc.i15no3<=0);
bgc.d15no3(idx)=nan;

idx = find(bgc.no2 < bgc.IsoThreshold | bgc.i15no2<=0);
bgc.d15no2(idx)=nan;

idx = find(bgc.nh4 < bgc.IsoThreshold | bgc.i15nh4<=0);
bgc.d15nh4(idx)=nan;

%ii = dNiso('i15N', bgc.i15no3, 'i14N', bgc.no3);
%idx = find(bgc.no3<bgc.IsoThreshold | bgc.i15no3<=0);
%bgc.d15no3 = ii.d15N;
%bgc.d15no3(idx)=nan;
%ii = dNiso('i15N', bgc.i15no2, 'i14N', bgc.no2);
%idx = find(bgc.no2<bgc.IsoThreshold | bgc.i15no2<=0);
%bgc.d15no2 = ii.d15N;
%bgc.d15no2(idx)=nan;
%ii = dNiso('i15N', bgc.i15nh4, 'i14N', bgc.nh4);
%idx = find(bgc.nh4<bgc.IsoThreshold | bgc.i15nh4<=0);
%bgc.d15nh4 = ii.d15N;
%bgc.d15nh4(idx)=nan;
%ii = dNiso('i15N', bgc.i15n2oA, 'i14N', bgc.n2o);
%idx = find(bgc.n2o<bgc.IsoThreshold/1000 | bgc.i15n2oA<=0);
%bgc.d15n2oA = ii.d15N;
%bgc.d15n2oA(idx)=nan;
%ii = dNiso('i15N', bgc.i15n2oB, 'i14N', bgc.n2o);
%idx = find(bgc.n2o<bgc.IsoThreshold/1000 | bgc.i15n2oB<=0);
%bgc.d15n2oB = ii.d15N;
%bgc.d15n2oB(idx)=nan;

end