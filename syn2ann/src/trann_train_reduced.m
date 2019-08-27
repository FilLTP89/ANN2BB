%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_train_: function to train ANN
%% *N.B.*
% Need for:
% _train_ann_justPSA.m_

%% *TRAIN ANN*
fprintf('---------------------\n1. TRAIN ANN\n---------------------\n');
for i_ = 1:ann.trn.nr
    train_ann_justPSA_reduced(ann.trn.wd,ann.trn.mtd(i_));
end

