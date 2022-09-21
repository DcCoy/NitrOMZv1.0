  function depth = bgc1d_detect_oxycline(o2,bgc,oxy_threshold);

 depth = nan(1,2);

 if nargin<3
    oxy_threshold = 1.0;
 end

 ind_anoxix = find(o2 < oxy_threshold);

 if length(ind_anoxix) >= 2
    depth(1) = bgc.zgrid(ind_anoxix(1));
    depth(2) = bgc.zgrid(ind_anoxix(end));
 elseif length(ind_anoxix) == 1
    depth(1) = bgc.zgrid(ind_anoxix(1));
    depth(2) = nan;
 elseif length(ind_anoxix) == 0
    depth(1) = nan;
    depth(2) = nan;
 end

