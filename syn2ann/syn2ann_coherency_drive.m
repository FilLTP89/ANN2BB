%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_coherency_drive_: function to compute cross-correlation between
% recorded and synthetic signals
%% *N.B.*
% Need for:
% __
fprintf('---------------------\n6. SPECTRAL MATCHING (JUST PSA)\n---------------------\n');
%% *MATCHING*
fprintf('--> Matching (just PSA)\n');
mon.cmb = nchoosek(1:mon.na,2)';
keyboard

switch lower(hybrid_type)
    case 'sp96'
        %
        % _SABETTA & PUGLIESE 1996_
        %
        spm.sps = syn2ann_sm(hbs.sps,trs.sps);
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        
        for j_ = 1:hbs.sps.mon.nc
            coh.rec.org.(hbs.sps.mon.cp{j_}) = ...
                syn2ann_coherency(mon,rec.org,(hbs.sps.mon.cp{j_}));
            coh.nss.org.(hbs.sps.mon.cp{j_}) = ...
                syn2ann_coherency(mon,nss.org,(hbs.sps.mon.cp{j_}));
            coh.hbs.sps.(hbs.sps.mon.cp{j_}) = ...
                syn2ann_coherency(mon,hbs.sps,(hbs.sps.mon.cp{j_}));
            coh.spm.sps.(hbs.sps.mon.cp{j_}) = ...
                syn2ann_coherency(mon,spm.sps.(hbs.sps.mon.cp{j_}),(hbs.sps.mon.cp{j_}));
        end
    case 'exsim'
        
    case 'both'
        
end
