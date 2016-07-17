%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. CORRECTION\n---------------------\n');
sps.org = [];
sps.hyb = [];
%% *RESAMPLING*
fprintf('--> Resampling\n');
nss.hyb = syn2ann_cornum(nss.org);
hbs     = nss.hyb;
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
nss.hyb = syn2ann_thp(nss.hyb);
hbs     = syn2ann_thp(hbs);
%% *SPECTRA*
fprintf('--> Spectra\n');
nss.hyb = syn2ann_spp(nss.hyb);
hbs = syn2ann_spp(hbs);