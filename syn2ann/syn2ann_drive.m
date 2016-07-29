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
% Need for
% _syn2ann_setup_emilia_num.m, syn2ann_records.m, syn2ann_numerical.m, 
% syn2ann_sp96.m, syn2ann_exsim.m, syn2ann_hybrid.m, syn2ann_justnum.m,
% syn2ann_ann.m, syn2ann_scaling.m, syn2ann_plot_res.m

%% *SET-UP*
syn2ann_setup_emilia;

%% *RECORDS*
syn2ann_records;

%% *NUMERICAL SIMULATIONS*
syn2ann_numerical;

if hybrid_flag
    switch lower(hybrid_type)
        case 'sp96'
            %% *SABETTA & PUGLIESE SYNTHETICS*
            syn2ann_sp96;
        case 'exsim'
            %% *SABETTA & PUGLIESE SYNTHETICS*
            syn2ann_exsim;
    end
    %% *LF-HF HYBRIDIZATION*
    syn2ann_hybrid;
else
    syn2ann_justnum;
end
%% *ANN - DATABASE*
syn2ann_ann;

%% SPECTRAL MATCHING
syn2ann_scaling;

%% PLOT RESULTS
syn2ann_plot_res;