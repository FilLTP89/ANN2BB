%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_scaling_drive_: function to perform spectral scaling upon ANN
% predictions
%% *N.B.*
% Need for:
% _syn2ann_sm.m,syn2ann_thp.m,syn2ann_spp.m_
fprintf('---------------------\n6. SPECTRAL MATCHING (JUST PSA)\n---------------------\n');
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