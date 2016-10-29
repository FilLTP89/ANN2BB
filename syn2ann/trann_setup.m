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

ccc;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%% *WORKDIR*
% _main workdir_
wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
    'EMILIA_2905');
sp = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
    'heavy_images');
ann.wd = fullfile(wd,'training');
fprintf('Workdir: %s\n',ann.wd);
% _database_
dbn = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
   'EMILIA_2905','simbad_v05','SIMBAD_v05_1.mat');
fprintf('Training Database: %s\n',dbn);

%% *LOAD ALL METADATA AVAILABLE*
trann_test_list;
% _select test case : selected_case = [a,b,...,d]_
selected_case = 1:4;

%% *DEFINE REAL RECORDS METADATA*
% _path to record files_
bhr.pt  = fullfile(wd,'records');
fprintf('--> Record Path: %s\n',bhr.pt);


bhr.ns = numel(selected_case);
for m_ = 1:bhr.ns
    for n_ = 1:fnn.bhrr
        bhr.(fni.bhrr{n_}){m_} = bhrr.(fni.bhrr{n_}){selected_case(m_)};
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

% _reference system_
bhr.rs = {'ew';'ns';'ud'};
bhr.cp = {'ew';'ns';'ud'};
bhr.nc = numel(bhr.cp);
[~,bhr.ci] = ismember(bhr.cp,bhr.rs);
fprintf('--> Reference system : ');
cellfun(@(x) fprintf('%s; ',x),bhr.rs);
fprintf('\n');
fprintf('--> Directions: ');
cellfun(@(x) fprintf('%s ',x),bhr.cp);
fprintf('\n');
% _motion components_
bhr.rc = {'a','v','d'};
bhr.nr = numel(bhr.rc);
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s; ',x),bhr.rc);
fprintf('\n');
bhr.cp = bhr.cp(:);
bhr.ci = bhr.ci(:);
fprintf('---------------------------------------------------------------\n');