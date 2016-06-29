%% *RECORDS*
fprintf('---------------------\n1. RECORDS\n---------------------\n');
%% *PARSING*
fprintf('--> Parsing\n');
% _original_
bhr.lfr = [];
bhr.hfr = [];
[bhr,rec.org]= rc_parser(bhr);
% % _filtered_
% bhr.lfr = [];
% bhr.hfr = 1.5;
% [~,rec.fil]= rc_parser(bhr);
%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
rec.org = syn2ann_thp(rec.org);
% rec.fil = syn2ann_thp(rec.fil);
%% *SPECTRA*
fprintf('--> Spectra\n');
rec.org = syn2ann_spp(rec.org);
% rec.fil = syn2ann_spp(rec.fil);