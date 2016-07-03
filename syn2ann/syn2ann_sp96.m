%% *SABETTA & PUGLIESE SYNTHETICS*
fprintf('---------------------\n3. SABETTA&PUGLIESE\n---------------------\n');
%% *METADATA*
mfn = fullfile(wd,'metadata.dat');
fprintf('--> Metadata: %s\n',mfn);
%% *GENERATION*
fprintf('--> Generation\n');
% _original_
mon.lfr = [];
mon.hfr = [];
sps.org = sp_generator(mon,mfn);
% % _filtered_
% mon.lfr = 1.5;
% mon.hfr = [];
% sps.fil = sp_generator(mon,mfn);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
sps.org = syn2ann_thp(sps.org);
% sps.fil = syn2ann_thp(sps.fil);
%% *SPECTRA*
fprintf('--> Spectra\n');
sps.org = syn2ann_spp(sps.org);
% sps.fil = syn2ann_spp(sps.fil);