close all; clc
clearvars -except hyb_EW_th hyb_NS_th hyb_Z_TH hyb_EW_RS hyb_NS_RS hyb_Z_RS 

%% Step 0: Defining input data
Setup;

%% Step 1: Loading the Input motions
load_the_motions;

%% Step 2: Calculating Synthetic with Sabetta and Pugliese '96                
Synthetic_SP96;

%% Step 3: LF-HF Hybridization             
LF_HF_Hybridization;

%% Step 4: Combining hybrid with ANN-DATABASE
ANN_Combination;

%% Step 5: Final Results 
Final_Results