%% *SET-UP*
% _general_
ccc;
fprintf('---------------------\n0. SETUP\n---------------------\n');
%%
% _workdir_
% wd = fullfile(filesep,'media','filippo','Data','Filippo','PHD_heavyweight',...
%     'EMILIA_2905');%fullfile(filesep,'media','filippo','Data','Filippo','PHD_passing_through','WORKDIR','HISADA','nco_eq','postprocess','hisada','shiba','nigh12');
wd = fullfile(filesep,'media','user','DATI','Filippo','PHD_heavyweight',...
    'EMILIA_2905');%fullfile(filesep,'media','user','DATI','Filippo','PHD_passing_through','WORKDIR','HISADA','nco_eq','postprocess','hisada','shiba','nigh12');
fprintf('Workdir: %s\n',wd);

%% *REAL RECORDS: METADATA*
% _path to record files_
bhr.pt  = fullfile(wd,'records');
fprintf('--> Record Path: %s\n',bhr.pt);
%%
% _station identity_
bhr.st{1}.id = {'MRN'};
% bhr.st{2}.id = {'MIR'};
bhr.st{1}.ni = {'IT'};
% bhr.st{2}.ni = {'TV'};
bhr.ns = numel(bhr.st);
fprintf('--> N. Stations: %u\n',bhr.ns);
cellfun(@(y) structfun(@(x) fprintf('--> Station ID: %s\n',x{:}),y,'UniformOutput',0),bhr.st);
fprintf('\n');
%%
% _recorded events_
bhr.ev = {'20120529.070003'};
bhr.ne  = numel(bhr.ev);
fprintf('--> N. Events: %u',bhr.ne);
fprintf('\n');
fprintf('--> Event IDs: ');
cellfun(@(x) fprintf('%s; ',x),bhr.ev);
fprintf('\n');
%%
% _database_
bhr.tp = {'itaca'};%{'kknpp'};
cellfun(@(x) fprintf('--> Database: %s',x),bhr.tp);
fprintf('\n');
%%
% _event label_
bhr.lb = {'Main Shock'};
%%
% _device list_
bhr.st{1}.dv = {''};
% bhr.st{2}.dv = {'01'};
for i_ = 1:bhr.ns
    bhr.nd(i_) = numel(bhr.st{i_}.dv);    
end

% _reference system_
bhr.rs = {'e';'n';'z'};
bhr.cp = {'n';'z'};
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

%% *MONITORS*
% _path to monitor files_
mon.pt  = fullfile(wd,'monitor');%fullfile(wd,'Copy_of_opm4');
fprintf('--> Monitor Path: %s\n',mon.pt);
% _metadata filename_
mon.fn  = fullfile(wd,'SM_Stations_Monitors.csv');%fullfile(wd,'SM_Stations_Monitors_Hisada.csv');
fprintf('--> Monitor File: %s\n',mon.fn);
% _type of simulation_
mon.tp  = 'S';%'H';
fprintf('--> Type of Simulation: %s\n',mon.tp);
% _monitor identity_
mon.id  = [16928];%[1,2];
mon.na  = numel(mon.id);
fprintf('--> N. Monitor: %u\n',mon.na);
arrayfun(@(x) fprintf('--> Monitor ID: %u; ',x),mon.id);
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
mon.rc  = {'d'};%{'v'};
mon.nr  = numel(mon.rc);
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s; ',x),mon.rc);
fprintf('\n');

bhr.cp = bhr.cp(:);
mon.cp = mon.cp(:);
bhr.ci = bhr.ci(:);
mon.ci = mon.ci(:);

%% *FILTERING RECORDS*
if strcmpi(mon.tp,'h')
    mon.fa = 1; % in Hz
    mon.fb = 1.5; % in Hz
else
    mon.fa = 1.3; % in Hz
    mon.fb = 1.7; % in Hz
end

%% *ANN*
ann.mtd.nl = {'net_075s_M.mat';'net_075s.mat';'net_1s.mat'};
ann.mtd.tc = {0.75;0.75;1};