%% *SPECTRAL SCALING*
fprintf('---------------------\n5. SPECTRAL MATCHING (JUST PSA)\n---------------------\n');
%% *MATCHING*
fprintf('--> Matching (just PSA)\n');
switch lower(hybrid_type)
    case 'sp96'
        %
        % _SABETTA & PUGLIESE 1996_
        %
        spm.sps = syn2ann_sm(hbs.sps,trs.sps);
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        
        for j_ = 1:hbs.sps.mon.nc
            spm.sps.(hbs.sps.mon.cp{j_}) = syn2ann_thp(spm.sps.(hbs.sps.mon.cp{j_}));
            spm.sps.(hbs.sps.mon.cp{j_}) = syn2ann_spp(spm.sps.(hbs.sps.mon.cp{j_}));
        end
    case 'exsim'
        %
        % _EXSIM_
        %
        spm.exs = syn2ann_sm(hbs.exs,trs.exs);
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        
        for j_ = 1:hbs.exs.mon.nc
            spm.exs.(hbs.exs.mon.cp{j_}) = syn2ann_thp(spm.exs.(hbs.exs.mon.cp{j_}));
            spm.exs.(hbs.exs.mon.cp{j_}) = syn2ann_spp(spm.exs.(hbs.exs.mon.cp{j_}));
        end
    case 'both'
        %
        % _SABETTA & PUGLIESE 1996_
        %
        spm.sps = syn2ann_sm(hbs.sps,trs.sps);
        %
        % _EXSIM_
        %
        spm.exs = syn2ann_sm(hbs.exs,trs.exs);
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        
        for j_ = 1:hbs.sps.mon.nc
            spm.sps.(hbs.sps.mon.cp{j_}) = syn2ann_thp(spm.sps.(hbs.sps.mon.cp{j_}));
            spm.sps.(hbs.sps.mon.cp{j_}) = syn2ann_spp(spm.sps.(hbs.sps.mon.cp{j_}));
        end
        
        for j_ = 1:hbs.exs.mon.nc
            spm.exs.(hbs.exs.mon.cp{j_}) = syn2ann_thp(spm.exs.(hbs.exs.mon.cp{j_}));
            spm.exs.(hbs.exs.mon.cp{j_}) = syn2ann_spp(spm.exs.(hbs.exs.mon.cp{j_}));
        end
end
