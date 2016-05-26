%%  *Spectral Matching: Numerical synthetics & ANN*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _synthetics2ann_matchup_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA) to target spectra obtained via
% Artificial Neural Networks.
%% INPUT:
%% OUTPUT:
%% N.B.
% Need for _ccc.m,plot_set_up.m,ns_parser.m,ns_spectra.m,sp_generator.m,
% sp_spectra.m,lfhf_mashup.m,nn_parser.m,ann2hbs_train.m,
% synthetics2ann_spectral_matching.m_

%% SET-UP
% _general_
ccc;[xtf,xlf,xtT,xlT]=plot_set_up;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%%
% _workdir_
wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
    'EMILIA_2905');
% wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_passing_through',...
%     'WORKDIR','HISADA','nco_eq','postprocess','hisada','shiba','nigh12');
fprintf('Workdir: %s\n',wd);
%% MONITORS
% _monitors_
% path to monitor files
mon.pt  = fullfile(wd,'monitor');
%mon.pt  = fullfile(wd,'Copy_of_opm4');
fprintf('--> Monitor Path: %s\n',mon.pt);
% metadata filename
mon.fn  = fullfile(wd,'SM_Stations_Monitors.csv');
%mon.fn  = fullfile(wd,'SM_Stations_Monitors_Hisada.csv');
fprintf('--> Monitor File: %s\n',mon.fn);
% type of monitor (S,H)
mon.tp  = 'S';
%mon.tp  = 'H';
fprintf('--> Type of Simulation: %s\n',mon.tp);
% monitor identity
mon.id  = [16928,3695];
%mon.id = [1,2];
% number of monitors
mon.na  = numel(mon.id);
fprintf('--> N. Monitor: %u\n',mon.na);
arrayfun(@(x) fprintf('--> Monitor ID: %u\n',x),mon.id);
% monitor recorded time-history (a,v,d)
mon.rc  = {'d'};
%mon.rc  = {'v'};
cellfun(@(x) fprintf('--> Components: %s\n',x),mon.rc);
% number of records
mon.nr  = numel(mon.rc);
% motion component (x,y,z)
mon.cp  = {'z','y','x'};
fprintf('--> Components: \n');
cellfun(@(x) fprintf('%s ',x),mon.cp);
fprintf('\n');
% number of components
mon.nc  = numel(mon.cp);

%% NUMERICAL SIMULATIONS
fprintf('---------------------\n1. NUMERICAL SIMULATIONS\n---------------------\n');
% _parse simulations_
fprintf('--> Parsing\n');
[mon,nss] = ns_parser(mon);
%%
% _compute spectra_
fprintf('--> Spectra\n');
nss = ns_spectra(nss);

%% SABETTA & PUGLIESE SYNTHETICS
fprintf('---------------------\n2. SABETTA&PUGLIESE\n---------------------\n');
% _generate SP synthetics_
fprintf('--> Generation\n');
sps = sp_generator(nss.mon,fullfile(wd,'metadata.dat'));
fprintf('--> Metadata: %s\n',sps.mon.fn);
%%
% _compute spectra_
fprintf('--> Spectra\n');
sps = sp_spectra(sps);

%% LF-HF HYBRIDIZATION
% _spectral mashup LF/HF_
hbs = lfhf_mashup(nss,sps);

%% ANN - DATABASE
% _ parse trained ANN network - horizontal directions_
% _parse trained ANN network - x direction_
if ismember({'x'},mon.cp)
    fn = fullfile(wd,'training','net_075s.mat');
    TnC  = 0.75;
    ann.x = nn_parser(TnC,fn);
end
% _parse trained ANN network - y direction_
if ismember({'y'},mon.cp)
    fn = fullfile(wd,'training','net_075s.mat');
    TnC  = 0.75;
    ann.y = nn_parser(TnC,fn);
end
% _parse trained ANN network - z direction_
if ismember({'z'},mon.cp)
    fn = fullfile(wd,'training','net_075s.mat');
    TnC  = 0.75;
    ann.z = nn_parser(TnC,fn);
end
ann.cp = fieldnames(ann);
[~,~,ib] = intersect(hbs.mon.cp,ann.cp,'stable');
ann.cp = ann.cp(ib);
% _ apply trained ANN on hybrid accelerograms_
trs = ann2hbs_train(hbs,ann);

%% SPECTRAL MATCHING
% _spectral scaling_ hybrid synthetics + ANN_
synthetics2ann_spectral_matching(hbs,trs);

%% POST-PROCESS
%%
% _hybridization - compare PSA-spectra_
fpplot('xpl',{nss.mon.vTn,sps.mon.vTn,hbs.mon.vTn},...
    'ypl',{nss.syn{1}.psa.x,sps.syn{1}.psa.x,hbs.syn{1}.psa.x,...
    nss.syn{1}.pga.x,sps.syn{1}.pga.x,hbs.syn{1}.pga.x},...
    'pax',{1;2;3},'spg',[1,3],'pfg',[0 0 20 8],...
    'scl',{'log','log','log'},'xlm',{[.1,10]},'xlb',{'T [s]'},...
    'ylm',{[1e-1,1e1]},'ylb',{'PSA [m/s/s]','',''},'tit',{'LF';'HF';'BB'});
%%
% _hybridization - single PSA-spectra_
fpplot('xpl',{nss.mon.vTn,sps.mon.vTn,hbs.mon.vTn},...
    'ypl',{abs(nss.syn{1}.psa.x),abs(sps.syn{1}.psa.x),...
    abs(hbs.syn{1}.psa.x)},'pfg',[0 0 16 16],'scl',{'log','log','log'},...
    'xlm',{[.1,10]},'xlb',{'T [s]'},'ylm',{[1e-1,1e1]},'ylb',{'PSA [m/s/s]','',''},...
    'leg',{'LF';'HF';'BB'},'tit',{'BB-PSA'});
%%
% _hybridization - compare FS-spectra_
fpplot('xpl',{nss.mon.vfr{1},sps.mon.vfr{1},hbs.mon.vfr{1}},...
    'ypl',{abs(nss.syn{1}.fsa.x),abs(sps.syn{1}.fsa.x),...
    abs(hbs.syn{1}.fsa.x)},'pax',{1;2;3},'spg',[1,3],'pfg',[0 0 20 8],...
    'scl',{'log','log','log'},'xlm',{[.1,10]},'xlb',{'f [Hz]'},...
    'ylm',{[1e-2,1e1]},'ylb',{'FS [1]','',''},'tit',{'LF';'HF';'BB'});
%%
% _hybridization - single FS-spectra_
fpplot('xpl',{nss.mon.vfr{1},sps.mon.vfr{1},hbs.mon.vfr{1}},...
    'ypl',{abs(nss.syn{1}.fsa.x),abs(sps.syn{1}.fsa.x),...
    abs(hbs.syn{1}.fsa.x)},'pfg',[0 0 16 16],'scl',{'log','log','log'},...
    'xlm',{[.1,10]},'xlb',{'f [Hz]'},'ylm',{[1e-2,1e1]},'ylb',{'FS [1]','',''},...
    'leg',{'LF';'HF';'BB'},'tit',{'BB-Fourier Spectra'});
%%
% _hybridization - compare time-histories_
fpplot('xpl',{nss.mon.vtm{1},sps.mon.vtm{1},hbs.mon.vtm{1}},...
    'ypl',{nss.syn{1}.tha.x,sps.syn{1}.tha.x,hbs.syn{1}.tha.x},...
    'pax',{1;2;3},'spg',[1,3],'pfg',[0 0 20 8],...
    'xlb',{'t [s]'},'ylb',{'a(t) [m/s/s]','',''},'tit',{'LF';'HF';'BB'});
%%
% _hybridization - single spectra_
fpplot('xpl',{nss.mon.vtm{1},sps.mon.vtm{1},hbs.mon.vtm{1}},...
    'ypl',{nss.syn{1}.tha.x,sps.syn{1}.tha.x,hbs.syn{1}.tha.x},...
    'pfg',[0 0 16 16],'xlb',{'t [s]'},'ylb',{'a(t) [m/s/s]','',''},...
    'leg',{'LF';'HF';'BB'},'tit',{'BB-Fourier Spectra'});



