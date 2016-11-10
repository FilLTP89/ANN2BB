%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_drive_: function to train and test ANN
%% *N.B.*
% Need for:
% _ccc.m,plot_set_up.m,ns_parser.m,ns_spectra.m,sp_generator.m,
% sp_spectra.m,lfhf_mashup.m,nn_parser.m,ann2hbs_train.m,
% synthetics2ann_spectral_matching.m_
%% *REFERENCES*
% @Book{Book_Haykin_1999_ANN,
%   Title                    = {{Neural Newtworks: a Comprehensive Foundation}},
%   Author                   = {Haykin, S.},
%   Publisher                = {Prentice-Hall International, Inc.},
%   Year                     = {1999},
% 
%   File                     = {Book_Haykin_1999_ANN.pdf:Book_Haykin_1999_ANN.pdf:PDF}
% }

%% *TRAIN SET-UP*
trann_setup;

%% *TRAIN ANN*
%trann_train;

%% *TEST TRAINED ANN*
trann_load;
trann_test;

%% *PLOTTING TEST RESULTS*
trann_test_plot_single;
%trann_test_plot_compare;
