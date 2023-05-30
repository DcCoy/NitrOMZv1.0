o2vals = linspace(33,200,10);

for i = 1:length(o2vals)
    o2 = o2vals(i);
    savestr = insertAfter("../saveout/o2sensitivity.mat","sensitivity",string(o2));
    disp(savestr)
    bgc = bgc1d_run('ParNames',"o2_bot", 'ParVal', o2);
    save(savestr,"bgc");
end