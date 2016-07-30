%% *SPECTRAL SCALING*
fprintf('---------------------\n5. SPECTRAL MATCHING (JUST PSA)\n---------------------\n');
%% *MATCHING*
fprintf('--> Matching (just PSA)\n');
spm = syn2ann_sm(hbs,trs);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
for j_ = 1:hbs.mon.nc
    spm.(hbs.mon.cp{j_}) = syn2ann_thp(spm.(hbs.mon.cp{j_}));
    spm.(hbs.mon.cp{j_}) = syn2ann_spp(spm.(hbs.mon.cp{j_}),2);
end