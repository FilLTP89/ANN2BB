%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_run_: function to perform classical LF/HF hybridization between
% SPEED results (PBS) and empirical/stochastic synthetics (SPS/EXS) and spectral
% match them upon ANN predictions.
%% *N.B.*
% Working-Flow:
% 1). PARSING PBS: parse physics-based simulation outputs
% 2). ANN SIMULATION: apply ANN to PBS
% 3). HYBRIDIZE WITH GOF: classical LF/HF hybridization + best PSA-GoF
% 4). SPECTRAL MATCHING: spectral matching best hybrid upon ANN predictions
% Need for:
% _syn2ann_pbs_drive.m, syn2ann_ann_drive.m, syn2ann_super_hybridize.m,
% syn2ann_scaling_drive.m_
% 
% #checked with Ali/Maria
%% *1). PARSING PBS*
syn2ann_pbs_drive;

% #checked with Ali/Maria
%% *2). ANN SIMULATION*
syn2ann_ann_drive;

%% *3). HYBRIDIZE WITH GOF*
syn2ann_super_hybridize;

%% *4). SPECTRAL MATCHING*
syn2ann_scaling_drive;
