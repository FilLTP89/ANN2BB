%% *RECORDS*
fprintf('---------------------\n1. RECORDS\n---------------------\n');
%% *PARSING*
fprintf('--> Parsing\n');
% _original_
ew = importdata(fullfile(bhr.pt,'EW_corrette.csv'));
ns = importdata(fullfile(bhr.pt,'NS_corrette.csv'));
ud = importdata(fullfile(bhr.pt,'UD_corrette.csv'));
cor.MRN.tha.e = ew.data(:,2);
cor.MRN.thv.e = ew.data(:,14);
cor.MRN.thd.e = ew.data(:,26);
cor.MRN.tha.n = ns.data(:,2);
cor.MRN.thv.n = ns.data(:,14);
cor.MRN.thd.n = ns.data(:,26);
cor.MRN.tha.z = ud.data(:,2);
cor.MRN.thv.z = ud.data(:,14);
cor.MRN.thd.z = ud.data(:,26);

cor.MIR08.tha.e = ew.data(:,6);
cor.MIR08.thv.e = ew.data(:,18);
cor.MIR08.thd.e = ew.data(:,30);
cor.MIR08.tha.n = ns.data(:,6);
cor.MIR08.thv.n = ns.data(:,18);
cor.MIR08.thd.n = ns.data(:,30);
cor.MIR08.tha.z = ud.data(:,6);
cor.MIR08.thv.z = ud.data(:,18);
cor.MIR08.thd.z = ud.data(:,30);
bhr.lfr = [];
bhr.hfr = [];
[bhr,rec.org]= syn2ann_rec_parser(bhr,cor);

%% *BASELINE CORRECTED*
rec.fil = syn2ann_blc(rec.org);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
rec.org = syn2ann_thp(rec.org);
rec.fil = syn2ann_thp(rec.fil);
%% *SPECTRA*
fprintf('--> Spectra\n');
rec.org = syn2ann_spp(rec.org);
rec.fil = syn2ann_spp(rec.fil);