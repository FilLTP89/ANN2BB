%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_rec_drive_: function to parse and correct recordings
%% *N.B.*
% Need for:
% _importdata.m,syn2ann_rec_parser.m,syn2ann_blc.m,syn2ann_thp.m,syn2ann_spp.m_
fprintf('---------------------\n1. RECORDS\n---------------------\n');
%% *PARSING RECORDS*
fprintf('--> Parsing\n');
% _original_
ew = importdata(fullfile(bhr.pt,'EW_corrette.csv'));
ns = importdata(fullfile(bhr.pt,'NS_corrette.csv'));
ud = importdata(fullfile(bhr.pt,'UD_corrette.csv'));
cor.MRN.tha.ew = ew.data(:,2);
cor.MRN.thv.ew = ew.data(:,14);
cor.MRN.thd.ew = ew.data(:,26);
cor.MRN.tha.ns = ns.data(:,2);
cor.MRN.thv.ns = ns.data(:,14);
cor.MRN.thd.ns = ns.data(:,26);
cor.MRN.tha.ud = ud.data(:,2);
cor.MRN.thv.ud = ud.data(:,14);
cor.MRN.thd.ud = ud.data(:,26);

cor.MIR08.tha.ew = ew.data(:,6);
cor.MIR08.thv.ew = ew.data(:,18);
cor.MIR08.thd.ew = ew.data(:,30);
cor.MIR08.tha.ns = ns.data(:,6);
cor.MIR08.thv.ns = ns.data(:,18);
cor.MIR08.thd.ns = ns.data(:,30);
cor.MIR08.tha.ud = ud.data(:,6);
cor.MIR08.thv.ud = ud.data(:,18);
cor.MIR08.thd.ud = ud.data(:,30);
bhr.lfr = [];
bhr.hfr = [];
[bhr,rec.org]= syn2ann_rec_parser(bhr,cor);

%% *BASELINE CORRECTION*
rec.fil = syn2ann_blc(rec.org);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
rec.org = syn2ann_thp(rec.org);
rec.fil = syn2ann_thp(rec.fil);
%% *SPECTRA*
fprintf('--> Spectra\n');
rec.org = syn2ann_spp(rec.org);
rec.fil = syn2ann_spp(rec.fil);