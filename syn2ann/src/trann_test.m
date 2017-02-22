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

%% *TEST ANN*
fprintf('---------------------\n2. TEST ANN\n---------------------\n');
for i_ = 1:tst.mtd.nr
    
    ann.tst{i_} = test_ann_justPSA(ann.tst{i_},rec.org);
    
end

