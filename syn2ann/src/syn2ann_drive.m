%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_drive_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA/SEM3D) to target spectra obtained via
% Artificial Neural Networks.
%% *N.B.*
% By setting flag_map=1, the user can run a fast analysis to plot single station
% time-histories/spectra (flag_map=0, demonstrative case) or run an analysis to
% obtain shake maps (flag_map=1, parallel implementation on remote cluster).
% Need for:
% _syn2ann_setup_maps.m,syn2ann_run_maps.m,
% syn2ann_setup_fast.m,syn2ann_run_fast.m_
ccc;

flag_map = 0; % flag to produce map data
flag_plot_results=0; % flag to plot results
flag_sensitivity=0;

if flag_map % write map data
    %% *1). CUSTOMIZE ANALYSIS SET-UP*
    syn2ann_setup_maps;
    
    %% *2). RUN SYN2ANN TO GET SHAKE MAPS (DNC)*
    syn2ann_run_maps;
else
    if flag_sensitivity
        
        %% *1). CUSTOMIZE ANALYSIS SET-UP*
        syn2ann_setup_fast;
        
        syn2ann_test_train;
        
        %% *2). PARSING REC (DNC)*
        syn2ann_rec_drive;
        
        %% *3). RUN SYN2ANN ON SINGLE STATIONS (DNC)*
        syn2ann_run_sensitivity;
        
        %% *5). SAVE RESULTS (DNC)*
        syn2ann_save_res;
        
        %% *6). PLOT RESULTS*
        if flag_plot_results
            syn2ann_plot_res_single;
        end
        
    else
        
        %% *1). CUSTOMIZE ANALYSIS SET-UP*
        syn2ann_setup_fast;
        
        %% *52). PARSING REC (DNC)*
        syn2ann_rec_drive;
        
        %% *3). RUN SYN2ANN ON SINGLE STATIONS (DNC)*
        syn2ann_run;
        
        %% *5). SAVE RESULTS (DNC)*
        syn2ann_save_res;
        
        %% *6). PLOT RESULTS*
        if flag_plot_results
            syn2ann_plot_res_single;
        end
        
    end
end
