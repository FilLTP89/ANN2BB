%% *SABETTA & PUGLIESE SYNTHETICS*
fprintf('---> SABETTA&PUGLIESE\n\n');

%% *GENERATION*
fprintf('--> Generation\n');
% _original_
mon.lfr = [];
mon.hfr = [];
sps.org{NIT} = syn2ann_sp96_generator(mon,mtd.sp96);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
sps.org{NIT} = syn2ann_thp(sps.org{NIT});

%% *SPECTRA*
fprintf('--> Spectra\n');
sps.org{NIT} = syn2ann_spp(sps.org{NIT});