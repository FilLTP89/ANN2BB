clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change log
% AGO on 18/12/2016: Loop is added to create and store n hybrid motions
% From those motions, the best candidates for E,W, and Z are selected.

% note: further cleaning of parameters may be needed, since major
% re-organization of the code together with insertion/creation of selection 
% routines are provided.

% current version of the code is v3.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 0: Defining input data
Setup;

%% Step 1: Loading the Input motions
load_the_motions;

for x=1:1
    clearvars -except wd wd_data wd_results cfr_record Mw scc R_epi ...
    inp_Tn tar_Tn record num_sim net_h net_v ...
    x hyb_EW_acc hyb_NS_acc hyb_Z_acc hyb_EW_PSA hyb_NS_PSA hyb_Z_PSA

%% Step 2: Calculating Synthetic with Sabetta and Pugliese '96                
Synthetic_SP96;

%% Step 3: LF-HF Hybridization             
LF_HF_Hybridization;

%   target_spectra;
%   if x==1
%     setup_score;
%   end
%   score;  
  
end

% %% Step 4: Combining hybrid with ANN-DATABASE
% ANN_Combination;
% 
% %% Step 5: Final Results 
% Final_Results_2format

