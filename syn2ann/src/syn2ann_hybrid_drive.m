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
[pbs.sps{NIT},sps.org{NIT}] = lfhf_rsmpl(pbs.org,sps.org{NIT});
pbs.sps{NIT}.mon = pbs.org.mon;

%% *PAD LF/HF*
fprintf('--> Padding/Tapering\n');
[pbs.sps{NIT},sps.org{NIT}] = lfhf_pad(pbs.sps{NIT},sps.org{NIT});

%% *ALIGN LF/HF*
fprintf('--> Align records\n');
[pbs.sps{NIT},sps.org{NIT}] = lfhf_shift(pbs.sps{NIT},sps.org{NIT});
pbs.sps{NIT} = syn2ann_thp(pbs.sps{NIT});
sps.org{NIT} = syn2ann_thp(sps.org{NIT});
pbs.sps{NIT} = syn2ann_spp(pbs.sps{NIT});
sps.org{NIT} = syn2ann_spp(sps.org{NIT});

%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Hybridization\n');
[pbs.hyb.sps{NIT},sps.hyb{NIT},hbs.sps{NIT}] = lfhf_mashup(pbs.sps{NIT},sps.org{NIT});

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
pbs.hyb.sps{NIT} = syn2ann_thp(pbs.hyb.sps{NIT});
sps.hyb{NIT} = syn2ann_thp(sps.hyb{NIT});
hbs.sps{NIT} = syn2ann_thp(hbs.sps{NIT});

%% *SPECTRA*
fprintf('--> Spectra\n');
pbs.hyb.sps{NIT} = syn2ann_spp(pbs.hyb.sps{NIT});
sps.hyb{NIT} = syn2ann_spp(sps.hyb{NIT});
hbs.sps{NIT} = syn2ann_spp(hbs.sps{NIT});
