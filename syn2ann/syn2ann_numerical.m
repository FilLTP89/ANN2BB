%% *NUMERICAL SIMULATIONS*
fprintf('---------------------\n2. NUMERICAL SIMULATIONS\n---------------------\n');
%% *PARSING*
fprintf('--> Parsing\n');
% _original_
mon.lfr = [];
mon.hfr = [];
[mon,nss.org]= ns_parser(mon);
% % _filtered_
% mon.lfr = [];
% mon.hfr = 1.5;
% [~,nss.fil]= ns_parser(mon);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
nss.org = syn2ann_thp(nss.org);
% nss.fil = syn2ann_thp(nss.fil);
%% *SPECTRA*
fprintf('--> Spectra\n');
[nss.org] = syn2ann_spp(nss.org);
% [nss.fil] = syn2ann_spp(nss.fil);