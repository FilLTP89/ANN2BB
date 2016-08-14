%% *SPECTRAL SCALING*
fprintf('---------------------\n5. SPECTRAL MATCHING (JUST PSA)\n---------------------\n');
%% *MATCHING*
fprintf('--> Matching (just PSA)\n');
spm = syn2ann_sm(hbs,trs);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
for j_ = 1:hbs.mon.nc
    spm.(hbs.mon.cp{j_}) = syn2ann_thp(spm.(hbs.mon.cp{j_}));
    spm.(hbs.mon.cp{j_}) = syn2ann_spp(spm.(hbs.mon.cp{j_}));
    if strcmpi(hbs.mon.cp{j_},'z')
        
        hold all;
        plot(trs.z.mon.vTn,trs.z.syn{1}.psa.z*100,'gd-');
        plot(spm.z.mon.vTn,spm.z.syn{1}.psa.z*100,'c+-');
        keyboard
    end
    
end