%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_pbs_drive_: function to parse PBS records
%% *N.B.*
% Need for:
% _syn2ann_pbs_parser.m,syn2ann_thp.m,syn2ann_spp.m_

fprintf('============================\n');
fprintf('------2. PB SIMULATIONS-----\n');
fprintf('============================\n');

%% *PARSING ORIGINAL NUMERICAL SIMULATIONS*
fprintf('--> Parsing\n');
mon.lfr = 0.05;
mon.hfr = [];
[mon,pbs.org]= syn2ann_pbs_parser(mon,bhr);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
pbs.org = syn2ann_thp(pbs.org);

%% *SPECTRA*
fprintf('--> Spectra\n');
[pbs.org] = syn2ann_spp(pbs.org);