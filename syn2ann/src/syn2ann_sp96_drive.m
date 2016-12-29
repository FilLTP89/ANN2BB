%% *SABETTA & PUGLIESE SYNTHETICS*
fprintf('---> SABETTA&PUGLIESE\n\n');
%% *METADATA*
mfn = fullfile(wd,'metadata.dat');
fprintf('--> Metadata: %s\n',mfn);
%% *GENERATION*
fprintf('--> Generation\n');
% _original_
mon.lfr = [];
mon.hfr = [];
sps.org = syn2ann_sp96_generator(mon,mtd.sp96);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
sps.org = syn2ann_thp(sps.org);

%% *SPECTRA*
fprintf('--> Spectra\n');
sps.org = syn2ann_spp(sps.org);