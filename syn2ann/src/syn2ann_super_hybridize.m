%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_super_hybridize_: function to hybridize low-frequency synthetics (PBS)
% with HF synthetics from Empirical (SPS) or Stochastic (EXS) methods. SPS is obtained
% via Sabetta-Pugliese-1996 approach, whereas EXS with EXSIM code. The best hybrid is 
% iteratively selected by computing a GoF on the PSA spectrum.
% *N.B.*
% Need for:
% _syn2ann_emp_stp_drive.m,syn2ann_hybrid_drive.m,syn2ann_gof_setup.m,syn2ann_gof_compute.m,
% syn2ann_gof_best.m_
%% *REFERENCES*
%@Article{Motazedian_Atkinson_2005,
%  Title                    = {{Stochastic finite-fault modeling based on a dynamic corner frequency}},
%  Author                   = {Motazedian, D. and Atkinson, G. M.},
%  Journal                  = {Bulletin of the Seismological Society of America},
%  Year                     = {2005},
%  Month                    = {June},
%  Number                   = {3},
%  Pages                    = {995-1010},
%  Volume                   = {95},
%}
%@Article{Sabetta_Pugliese_1996,
%  Title                    = {{Estimation of Response Spectra and Simulation of Nonstationary Earthquake Ground Motions}},
%  Author                   = {Sabetta, F. and Pugliese, A.},
%  Journal                  = {Bulletin of the Seismological Society of America},
%  Year                     = {1996},
%  Number                   = {2},
%  Pages                    = {337-352},
%  Volume                   = {86},
%}

for NIT=1:MAXIT
    
    %% *GENERATE EMPIRICAL - PARSE STOCHASTIC*
    syn2ann_emp_sto_drive;
    
    %% *LF-HF CLASSIC HYBRIDIZATION*
    syn2ann_hybrid_drive;
    
    %% *SCORE THE HYBRID TIME HISTORIES*
    if NIT==1
        syn2ann_gof_setup;
    end
    syn2ann_gof_compute;
    
end

%% *COMPUTE BEST GOF*
syn2ann_gof_best;
