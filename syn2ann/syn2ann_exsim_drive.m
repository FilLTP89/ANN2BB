%% *SABETTA & PUGLIESE SYNTHETICS*
fprintf('---------------------\n4. EXSIM\n---------------------\n');
%% *METADATA*
mfn = fullfile(wd,'metadata.dat');
fprintf('--> Metadata: %s\n',mfn);
%% *GENERATION*
fprintf('--> Generation\n');
% _original_
mon.lfr = [];
mon.hfr = [];
exs.org = syn2ann_exsim_parser(mon,mtd);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
exs.org = syn2ann_thp(exs.org);

%% *SPECTRA*
fprintf('--> Spectra\n');
exs.org = syn2ann_spp(exs.org);