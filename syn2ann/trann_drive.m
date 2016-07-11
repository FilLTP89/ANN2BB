%%  *Spectral Matching: Numerical synthetics & ANN*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _trann_drive_: function to train ANN
%% N.B.
% Need for _ccc.m,plot_set_up.m,ns_parser.m,ns_spectra.m,sp_generator.m,
% sp_spectra.m,lfhf_mashup.m,nn_parser.m,ann2hbs_train.m,
% synthetics2ann_spectral_matching.m_

%% *SET-UP*
trann_setup;

%% *TRAIN-ANN*
trann_train;