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
fprintf('============================\n');
fprintf('----9. SPECTRAL MATCHING----\n');
fprintf('============================\n');

%% *SPECTRAL MATCHING*
spm.sps = syn2ann_sm(hbs.bst,trs.sps);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values/Arias/Spectra\n');
for j_ = 1:hbs.bst.mon.nc
    spm.sps.(hbs.bst.mon.cp{j_}) = syn2ann_thp(spm.sps.(hbs.bst.mon.cp{j_}));
    spm.sps.(hbs.bst.mon.cp{j_}) = syn2ann_spp(spm.sps.(hbs.bst.mon.cp{j_}));
end