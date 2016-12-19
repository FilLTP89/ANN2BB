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
fprintf('---------------------\n3. EMPIRICAL-STOCHASTIC\n---------------------\n');
syn2ann_sp96_drive;

% switch lower(hybrid_type)
%     case 'sp96'
%         % _SABETTA & PUGLIESE 1996_
%     case 'exsim'
%         % _EXSIM_
%         syn2ann_exsim_drive;
%     case 'both'
%         % _SABETTA & PUGLIESE 1996_
%         syn2ann_sp96_drive;
%         % _EXSIM_
%         syn2ann_exsim_drive;
% end