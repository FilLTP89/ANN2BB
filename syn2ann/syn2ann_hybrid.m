%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');
%% *RESAMPLING*
[nss.org,sps.org] = lfhf_rsmpl(nss.org,sps.org);
%% *PAD LF/HF*
[nss.org,sps.org] = lfhf_pad(nss.org,sps.org);
%% *TUKEY LF/HF*
[nss.org,sps.org] = lfhf_taper(nss.org,sps.org);
%% *ALIGN LF/HF*
[nss.org,sps.org] = lfhf_shift(nss.org,sps.org);
nss.org = syn2ann_thp(nss.org);
sps.org = syn2ann_thp(sps.org);
nss.org = syn2ann_spp(nss.org);
sps.org = syn2ann_spp(sps.org);
%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Generation\n');
[nss.hyb,sps.hyb,hbs] = lfhf_mashup(nss.org,sps.org);
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