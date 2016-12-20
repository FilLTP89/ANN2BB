%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_hybrid_drive_: function to perform classical LF/HF hybridization.
%% *N.B.*
% Need for:
% _lfhf_rsmpl.m,lfhf_pad.m,lfhf_shift.m,syn2ann_thp.m,syn2ann_spp.m,
% lfhf_mashup.m_
%% *LF-HF HYBRIDIZATION*
fprintf('============================\n');
fprintf('--------5. HYBRIDIZE--------\n');
fprintf('============================\n');

%% *RESAMPLING*
fprintf('--> Resampling\n');
[pbs.sps,sps.org] = lfhf_rsmpl(pbs.org,sps.org);

%% *PAD LF/HF*
fprintf('--> Padding/Tapering\n');
[pbs.sps,sps.org] = lfhf_pad(pbs.sps,sps.org);

%% *ALIGN LF/HF*
fprintf('--> Align records\n');
[pbs.sps,sps.org] = lfhf_shift(pbs.sps,sps.org);
pbs.sps = syn2ann_thp(pbs.sps);
sps.org = syn2ann_thp(sps.org);
pbs.sps = syn2ann_spp(pbs.sps);
sps.org = syn2ann_spp(sps.org);

%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Hybridization\n');
[pbs.hyb.sps,sps.hyb,hbs.sps] = lfhf_mashup(pbs.sps,sps.org);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
pbs.hyb.sps = syn2ann_thp(pbs.hyb.sps);
sps.hyb = syn2ann_thp(sps.hyb);
hbs.sps = syn2ann_thp(hbs.sps);

%% *SPECTRA*
fprintf('--> Spectra\n');
pbs.hyb.sps = syn2ann_spp(pbs.hyb.sps);
sps.hyb = syn2ann_spp(sps.hyb);
hbs.sps = syn2ann_spp(hbs.sps);