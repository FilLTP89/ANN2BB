%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_setup_: function to set-up the ANN training 
%% *N.B.*
% Need for:
% _ccc.m,_

fprintf('---------------------\n0. SETUP\n---------------------\n');
%% *WORKDIR*
% _main workdir_
wd = fullfile(filesep,'home','filippo','Data','Filippo',...
    'PHD_passing_through_polimi','syn2ann','database');
% _save path_
sp = fullfile(filesep,'home','filippo','Data','Filippo','aeolus',...
    'image_database1');
% _database path_
dbn = fullfile(filesep,'home','filippo','Data','Filippo',...
    'PHD_passing_through_heavyweight','simbad_v05','SIMBAD_v05_1_vigne.mat');

if exist(wd,'dir')~=7
    wd  = '/mssmat2/home/gattif/Documents/PHD_passing_through_polimi/syn2ann/database';
    sp  = '/tmp1/gattif/ann_train';
    dbn = '/mssmat2/home/gattif/Documents/SIMBAD_v05_1.mat';  
    spp = '/tmp1/gattif/ann_train';
end

ann.trn.wd = fullfile(wd,'training');
fprintf('Training Workdir: %s\n',ann.trn.wd);
fprintf('Training Database: %s\n',dbn);

%wkd     = '~/Documents/ares/workdir/ANN2BB/sensitivity';
%mat     = load(fullfile(wkd,'data_rho_fast_1.mat'));
%res.fnm = load(fullfile(wkd,'res_rho_fast_1.mat'));

trann_test_list_sensitivity_semNsobol;
run_selcase = 1:12;

%% *DEFINE REAL RECORDS METADATA*
% _path to record files_
bhr.pt  = fullfile(wd,'records');
fprintf('--> Record Path: %s\n',bhr.pt);

bhr.ns = numel(run_selcase);
for m_ = 1:bhr.ns
    for n_ = 1:fnn.bhrr
        bhr.(fni.bhrr{n_}){m_} = bhrr.(fni.bhrr{n_}){run_selcase(m_)};
    end
end

fprintf('--> N. Stations: %u\n',bhr.ns);
cellfun(@(y) structfun(@(x) fprintf('--> Station ID: %s\n',x{:}),y,'UniformOutput',0),bhr.st);
fprintf('---------------------------------------------------------------\n');

bhr.na = 0;
for i_ = 1:bhr.ns
    bhr.nd(i_) = numel(bhr.st{i_}.dv);
    bhr.na = bhr.na+bhr.nd(i_);
end

%% *GENERAL SET-UP (CUSTOMIZE)*

% _reference system_
bhr.rs = {'ew';'ns';'ud'};
% _select component (ew,ns,ud,gh)
bhr.cp = {'ew';'ns'};
bhr.nc = numel(bhr.cp);
[~,bhr.ci] = ismember(bhr.cp,bhr.rs);
fprintf('--> Reference system : ');
cellfun(@(x) fprintf('%s; ',x),bhr.rs);
fprintf('\n');
fprintf('--> Directions: ');
cellfun(@(x) fprintf('%s ',x),bhr.cp);
fprintf('\n');
% _motion components_
bhr.rc = {'a'};
bhr.nr = numel(bhr.rc);
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s; ',x),bhr.rc);
fprintf('\n');
bhr.cp = bhr.cp(:);
bhr.ci = bhr.ci(:);
fprintf('---------------------------------------------------------------\n');


%% *DEFINE SPEED MONITORS METADATA (CUSTOMIZE)*
% _path to monitor files_
% % # paper BSSA
% mon.pt  = fullfile(wd,'monitors');
% # KKNPP AS8
% FOLDED
% mon.pt  = ['/home/filippo/Data/Filippo/aeolus/SEM_calculs/ncoeq2007_small_scale_OK/ncoeq2007_as8/',...
%     'kknpp_as8_topo_rf5_small_bsn_wdmp_npml_sb88_grd_tsuda_water_fold'];
% # KKNPP AS4
% FOLDED
mon.pt  = ['/tmp1/gattif/SEM_calculs/ncoeq2007_small_scale_OK/ncoeq2007_as4/',...
    'kknpp_as4_topo_rf5_small_bsn_wdmp_npml_sb35_25_grd_tsuda_water_fold'];
% LAYERED
% % % % mon.pt  = ['/home/filippo/Data/Filippo/aeolus/SEM_calculs/ncoeq2007_small_scale_OK/ncoeq2007_as4/',...
% % % %     'kknpp_as4_topo_rf5_small_bsn_wdmp_npml_sb35_25_grd_tsuda_water_1d'];

fprintf('--> Monitor Path: %s\n',mon.pt);
% _metadata filename (for Sabetta-Pugliese)_
% % # paper BSSA
% mon.fnm  = fullfile(wd,'SM_Stations_Monitors.csv');
% # KKNPP
mon.fnm  = fullfile(wd,'SM_Stations_Monitors_sem3d.csv');
fprintf('--> Monitor File: %s\n',mon.fnm);

%% *TYPE OF SIMULATION (CUSTOMIZE)*
% % SPEED
% mon.typ = 'speed';

% SEM3D
mon.typ  = 'sem3d';
fprintf('--> Type of Simulation: %s\n',mon.typ);

% _monitor identity_
mon.na = bhr.ns;
for m_ = 1:mon.na
    for n_ = 1:fnn.monn
        mon.(fni.monn{n_})(m_) = monn.(fni.monn{n_})(run_selcase(m_));
    end
end

fprintf('--> N. Monitor: %u\n',mon.na);
arrayfun(@(x) fprintf('--> Monitor ID: %u \n',x),mon.id);
% _reference system_
mon.rs = bhr.rs;
mon.cp = bhr.cp;
mon.ci = bhr.ci;
mon.nc = numel(mon.cp);
fprintf('--> Directions: ');
cellfun(@(x) fprintf('%s ',x),mon.cp);
% _motion components_
mon.rc  = {'d'};
mon.nr  = numel(mon.rc);
fprintf('\n');
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s ',x),mon.rc);
fprintf('\n---------------------------------------------------------------\n');
mon.cp = mon.cp(:);
mon.ci = mon.ci(:);
mon.hyb = 'butter';

%% *MAP DATA (CUSTOMIZE)*
mon.map.flg = 0;
mon.map.typ = 'PGA';
mon.map.stk = 95;
mon.map.vTn.psa = 0.5;
mon.map.vTn.rsd = 0.75;
mon.map.fnm = '/tmp1/gattif/syn2ann_maps/map_mrn';

ann.mtd.scl = {'ALL';'ALL';'ALL'}; 
% _corner period for each ANN_
ann.mtd.TnC = {0.75;0.75;0.75};
% _ANN motion component : gh,ud (geometric mean horizontal, updip)_
ann.mtd.cpn = {'gh';'gh';'ud'};
for i_ = 1:numel(ann.mtd.TnC)
    ann.mtd.nl(i_) = {sprintf('net_%u_%s_%s_30n.mat',...
        round(ann.mtd.TnC{i_}*100),ann.mtd.scl{i_},ann.mtd.cpn{i_})};
    ann.mtd.tol(i_).low.psa = 0.1;
    ann.mtd.tol(i_).low.pga = 0.1;
    ann.mtd.tol(i_).hgh.psa = 0.1;
    ann.mtd.tol(i_).hgh.pga = 0.1;
end
% _tolerances_
ann.mtd.nit = [20;20;50];



