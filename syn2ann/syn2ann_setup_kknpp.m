%% *SET-UP - NCO 2007*
ccc;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%% *WORKDIR*
wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_passing_through',...
    'WORKDIR','HISADA','nco_eq','postprocess','hisada','shiba','nigh12');
% wd = fullfile(filesep,'media','user','DATI','Filippo','PHD_passing_through',...
%     'WORKDIR','HISADA','nco_eq','postprocess','hisada','shiba','nigh12');
fprintf('Workdir: %s\n',wd);
%% *REAL RECORDS: METADATA*
% _path to record files_
bhr.pt  = fullfile(filesep,'media','filippo','Data','Filippo','PHD_passing_through',...
    'WORKDIR','HISADA','nco_eq','records');
% bhr.pt  = fullfile(filesep,'media','user','DATI','Filippo','PHD_passing_through',...
%     'WORKDIR','HISADA','nco_eq','records');
fprintf('--> Record Path: %s\n',bhr.pt);
% _station identity_
bhr.st{1}.id = {'KSH'};
bhr.ns = numel(bhr.st);
fprintf('--> N. Stations: %u\n',bhr.ns);
cellfun(@(y) structfun(@(x) fprintf('--> Station ID: %s\n',x{:}),y,'UniformOutput',0),bhr.st);
fprintf('\n');
% _recorded events_
bhr.ev = {'0707161013'};
bhr.ne  = numel(bhr.ev);
fprintf('--> N. Events: %u',bhr.ne);
fprintf('\n');
fprintf('--> Event IDs: ');
cellfun(@(x) fprintf('%s; ',x),bhr.ev);
fprintf('\n');
% _database_
bhr.tp = {'kknpp'};
cellfun(@(x) fprintf('--> Database: %s',x),bhr.tp);
fprintf('\n');
% _event label_
bhr.lb = {'Main Shock'};
% _device list_
bhr.st{1}.dv = {'SG4'};
bhr.na = 0;
for i_ = 1:bhr.ns
    bhr.nd(i_) = numel(bhr.st{i_}.dv);
    bhr.na = bhr.na+bhr.nd(i_);
end
% _reference system_
bhr.rs = {'e';'n';'z'};
bhr.cp = {'e';'n';'z'};
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
cellfun(@(x) fprintf('%s ',x),bhr.rc);
fprintf('\n');
bhr.cp = bhr.cp(:);
bhr.ci = bhr.ci(:);
%% *MONITORS*
% _path to monitor files_
mon.pt  = fullfile(wd,'Copy_of_opm4');
fprintf('--> Monitor Path: %s\n',mon.pt);
% _metadata filename_
mon.fn  = fullfile(wd,'SM_Stations_Monitors_Hisada.csv');
fprintf('--> Monitor File: %s\n',mon.fn);
% _type of simulation_
mon.tp  = 'H';
fprintf('--> Type of Simulation: %s\n',mon.tp);
% _monitor identity_
mon.id  = [1];
mon.na  = numel(mon.id);
fprintf('--> N. Monitor: %u\n',mon.na);
arrayfun(@(x) fprintf('--> Monitor ID: %u ',x),mon.id);
fprintf('\n');
% _reference system_
mon.rs = bhr.rs;
mon.cp = bhr.cp;
mon.ci = bhr.ci;
mon.nc = numel(mon.cp);
fprintf('--> Directions: ');
cellfun(@(x) fprintf('%s ',x),mon.cp);
fprintf('\n');
% _motion components_
mon.rc  = {'v'};
mon.nr  = numel(mon.rc);
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s; ',x),mon.rc);
fprintf('\n');
mon.cp = mon.cp(:);
mon.ci = mon.ci(:);
%% *FILTERING RECORDS*
mon.fa = 1; % in Hz
mon.fb = 1.5; % in Hz
%% *ANN*
ann.mtd.nl = {'net_05s.mat';'net_075s.mat';'net_1s.mat'};
ann.mtd.tc = {0.5;0.75;1};