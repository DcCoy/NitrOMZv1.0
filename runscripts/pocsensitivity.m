pocvals = linspace(-15,-6,10);

for i = 1:length(pocvals)
    poc = pocvals(i);
    savestr = insertAfter("../saveout/pocsensitivity.mat","sensitivity",string(poc));
    disp(savestr)
    bgc = bgc1d_run('ParNames',"poc_flux_top", 'ParVal', poc./86400);
    save(savestr,"bgc");
end