%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_drive_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA) to target spectra obtained via
% Artificial Neural Networks.
%% *N.B.*
% Need for:
% _syn2ann_setup.m,syn2ann_rec_drive.m, syn2ann_sim_drive.m,
% syn2ann_emp_sto_drive.m, syn2ann_hybrid_drive.m, syn2ann_ann_drive.m, 
% syn2ann_scaling_drive.m, syn2ann_save_res, syn2ann_plot_res_single.m_
%% *REFERENCES*

%% *ANALYSIS SET-UP*
syn2ann_setup;

%% *PARSE RECORDS*
syn2ann_rec_drive;

%% *PARSE NUMERICAL SIMULATIONS*
syn2ann_sim_drive;

%% *GENERATE EMPIRICAL - PARSE STOCHASTIC*
syn2ann_emp_sto_drive;

%% *LF-HF CLASSIC HYBRIDIZATION*
syn2ann_hybrid_drive;

%% *PARSE ANN DATABASE*
syn2ann_ann_drive;

%% *HYB-ANN SPECTRAL MATCHING*
syn2ann_scaling_drive;

%% *HYB-ANN SPECTRAL MATCHING*
% syn2ann_coherency_drive;

%% *SAVE RESULTS*
%syn2ann_save_res;

%% *PLOT RESULTS*
syn2ann_plot_res_single;