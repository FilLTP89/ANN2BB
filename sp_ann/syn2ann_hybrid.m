%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');
%% *RESAMPLING*
[nss.hyb,sps.hyb] = lfhf_rsmpl(nss.org,sps.org);
%% *PAD LF/HF*
[nss.hyb,sps.hyb] = lfhf_pad(nss.hyb,sps.hyb);
%% *ALIGN LF/HF*
[nss.hyb,sps.hyb] = lfhf_shift(nss.hyb,sps.hyb);
nss.hyb = syn2ann_spp(nss.hyb);
sps.hyb = syn2ann_spp(sps.hyb);
%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Generation\n');
[nss.hyb,sps.hyb,hbs] = lfhf_mashup(nss.hyb,sps.hyb);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
nss.hyb = syn2ann_thp(nss.hyb);
sps.hyb = syn2ann_thp(sps.hyb);
hbs = syn2ann_thp(hbs);
%% *SPECTRA*
fprintf('--> Spectra\n');
nss.hyb = syn2ann_spp(nss.hyb);
sps.hyb = syn2ann_spp(sps.hyb);
hbs = syn2ann_spp(hbs);