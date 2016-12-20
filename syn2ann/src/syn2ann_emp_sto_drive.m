%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_emp_sto_drive_: function to generate stochastic and empirical
% wave-forms for classical LF/HF hybridization.
%% *N.B.*
% Need for:
% _syn2ann_sp96_drive.m,syn2ann_exsim_drive.m_
fprintf('============================\n');
fprintf('---------4. EMP/STO---------\n');
fprintf('============================\n');

%% *GENERATE EMPIRICAL RECORDS*
fprintf('--> Generate SP96 \n');
syn2ann_sp96_drive;