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
%%
% _general_
ccc;[xtf,xlf,xtT,xlT]=plot_set_up;
%%
% _workdir_
% wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
%     'EMILIA_2905');
wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_passing_through',...
    'WORKDIR','HISADA','nco_eq','postprocess','hisada','shiba','nigh12');

%% PARSING RECORDS
%%
% _monitors_
% path to monitor files
%mon.pt  = fullfile(wd,'monitor');
mon.pt  = fullfile(wd,'Copy_of_opm4');
% metadata filename
%mon.fn  = fullfile(wd,'SM_Stations_Monitors.csv');
mon.fn  = fullfile(wd,'SM_Stations_Monitors_Hisada.csv');
% type of monitor (S,H)
% mon.tp  = 'S';
mon.tp  = 'H';
% monitor identity
% mon.id  = [3695;16928];
mon.id = [1,2];
% number of monitors
mon.na  = numel(mon.id);                       
% monitor recorded time-history (a,v,d)
% mon.rc  = {'d'};
mon.rc  = {'v'};
% number of records
mon.nr  = numel(mon.rc);                        
% motion component (x,y,z)
mon.cp  = {'x','y','z'};
% number of components
mon.nc  = numel(mon.cp);
%%
% _parse simulations_
[mon,nss] = ns_parser(mon);
%%
% _compute spectra_
nss = ns_spectra(nss);
%% SABETTA & PUGLIESE SYNTHETICS
%%
% _generate SP synthetics_
sps = sp_generator(nss.mon,fullfile(wd,'metadata.dat'));
%% 
% _compute spectra_
sps = sp_spectra(sps);
%% LF-HF HYBRIDIZATION
%% 
% _spectral mashup LF/HF_
hbs = lfhf_mashup(nss,sps);
%% ANN - DATABASE
%%
% _parse trained ann network_
ann = nn_parser(0.75,fullfile(wd,'training','net_075s.mat'));
%%
% _train ANN on hybrid accelerograms_
trs = ann2hbs_train(hbs,ann);
%% SPECTRAL MATCHING
%%
% _spectral scaling_ hybrid synthetics + ANN_
synthetics2ann_spectral_matching(hbs,trs);