%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_test_: function to test trained ANN
%% *N.B.*
% Need for:
% _trann_rec_parser.m,test_ann_justPSA.m_

%% *LOAD TEST SET*
fprintf('---------------------\n1. LOAD TEST SET\n---------------------\n');
[bhr,rec]= trann_rec_parser(bhr);
rec = syn2ann_thp(rec);
rec = syn2ann_spp(rec);

%% *TEST ANN*
fprintf('---------------------\n1. TEST ANN\n---------------------\n');
for i_ = 1:tst.mtd.nr
     ann.tst{j_} = test_ann_justPSA(ann.tst{j_},rec);
end