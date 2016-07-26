%% *SPECTRAL SCALING*
fprintf('---------------------\n5. SPECTRAL MATCHING (WITH PGV)\n---------------------\n');
%% *MATCHING*
fprintf('--> Matching (PGV included)\n');
spm_withPGV = syn2ann_sm(hbs,trs_withPGV);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
for j_ = 1:hbs.mon.nc
    spm_withPGV.(hbs.mon.cp{j_}) = syn2ann_thp(spm_withPGV.(hbs.mon.cp{j_}));
    spm_withPGV.(hbs.mon.cp{j_}) = syn2ann_spp(spm_withPGV.(hbs.mon.cp{j_}),2);
end

%% *MATCHING*
fprintf('--> Matching (no PGV included)\n');
spm_noPGV = syn2ann_sm(hbs,trs_noPGV);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
for j_ = 1:hbs.mon.nc
    spm_noPGV.(hbs.mon.cp{j_}) = syn2ann_thp(spm_noPGV.(hbs.mon.cp{j_}));
    spm_noPGV.(hbs.mon.cp{j_}) = syn2ann_spp(spm_noPGV.(hbs.mon.cp{j_}),2);
end