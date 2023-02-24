function Data = proc_data(bgc,data_raw,tracers)

% Processes data with a spline to cover the model grid
Data.name = tracers;
Data.val=nan(length(tracers),bgc.nz);

fieldname = fields(data_raw);
for i = 1 : length(tracers)
   if sum(strcmp(tracers{i},fieldname)) == 1
      Data.val(i,:) =  data_raw.(tracers{i});
   elseif sum(strcmp(tracers{i},fieldname)) == 0 
      Data.val(i,:) = nan(1,bgc.nz);
   else
      error('strcmp(tracers,fieldname) should be 1 or 0');
   end
end

