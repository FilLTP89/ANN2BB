%% *NUMERICAL SIMULATIONS*
fprintf('---------------------\n2. NUMERICAL SIMULATIONS\n---------------------\n');
%% *PARSING*
fprintf('--> Parsing\n');
% _original_
mon.lfr = 0.05;
mon.hfr = [];
[mon,nss.org]= ns_parser(mon,bhr);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
nss.org = syn2ann_thp(nss.org);
%% *SPECTRA*
fprintf('--> Spectra\n');
[nss.org] = syn2ann_spp(nss.org);