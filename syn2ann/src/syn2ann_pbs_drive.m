%% *NUMERICAL SIMULATIONS*
fprintf('---------------------\n2. NUMERICAL SIMULATIONS\n---------------------\n');
%% *PARSING*
fprintf('--> Parsing\n');
% _original_
mon.lfr = 0.05;
mon.hfr = [];
%% PARSE ORIGINAL PBS RECORDS
[mon,pbs.org]= syn2ann_pbs_parser(mon,bhr);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
pbs.org = syn2ann_thp(pbs.org);

%% *SPECTRA*
fprintf('--> Spectra\n');
[pbs.org] = syn2ann_spp(pbs.org);