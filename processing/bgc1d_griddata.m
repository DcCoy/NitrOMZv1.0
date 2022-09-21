function [GridData] = ToDensToVert2(data,zgrid,bgc)

 rres  = (bgc.ztop - bgc.zbottom)/bgc.npt;
 rstep = (bgc.ztop - bgc.zbottom)/rres;
 Z = linspace(bgc.ztop,bgc.zbottom,rstep+1);
 clear tmp
 tmp(bgc.npt+1).val = [];
 for k = 1 : length(zgrid)
    if (zgrid(k) >= (bgc.zbottom - 0.5 * rres)) & (zgrid(k) <= (bgc.ztop + 0.5 * rres))
       [~, idx] = min(abs(zgrid(k) - Z));
       if ~isnan(zgrid(k))            
          if isempty(tmp(idx).val)
             tmp(idx).val = data(k);                     
          else
             tmp(idx).val = [tmp(idx).val data(k)];
          end
       end
    end
 end

 for t = 1 : bgc.npt+1
    GridData(t) = nanmean(tmp(t).val);
 end

