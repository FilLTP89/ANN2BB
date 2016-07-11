%%  *Spectral Matching: Numerical synthetics & ANN*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_drive_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA) to target spectra obtained via
% Artificial Neural Networks.
%% N.B.
% Need for _ccc.m,plot_set_up.m,ns_parser.m,ns_spectra.m,sp_generator.m,
% sp_spectra.m,lfhf_mashup.m,nn_parser.m,ann2hbs_train.m,
% synthetics2ann_spectral_matching.m_

%% *SET-UP*
syn2ann_setup_emilia;
% syn2ann_setup_kknpp;

%% *RECORDS*
syn2ann_records;

%% *NUMERICAL SIMULATIONS*
syn2ann_numerical;

%% *SABETTA & PUGLIESE SYNTHETICS*
syn2ann_sp96;

%% *LF-HF HYBRIDIZATION*
syn2ann_hybrid;

%% *ANN - DATABASE*
syn2ann_ann;

%% SPECTRAL MATCHING
syn2ann_scaling;

%% PLOT RESULTS
cd(wd);
syn2ann_plot_res_MRN;
% syn2ann_plot_res_MIR08;
% syn2ann_plot_res_AQK;
% syn2ann_plot_res_AQU;