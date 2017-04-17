clear bhr mon mtd
%% *DEFINE REAL RECORDS METADATA*
% _path to record files_
bhr.pt  = fullfile(wd,'records');
fprintf('--> Record Path: %s\n',bhr.pt);

bhr.ns = numel(selected_case);
ns = bhr.ns;
for m_ = 1:ns
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
bhr.rc = {'a','v','d'};
bhr.nr = numel(bhr.rc);
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s; ',x),bhr.rc);
fprintf('\n');
bhr.cp = bhr.cp(:);
bhr.ci = bhr.ci(:);
fprintf('---------------------------------------------------------------\n');

%% *DEFINE SPEED MONITORS METADATA*
% _path to monitor files_
mon.pt  = fullfile(wd,'monitors');
fprintf('--> Monitor Path: %s\n',mon.pt);
% _metadata filename_
mon.fnm  = fullfile(wd,'SM_Stations_Monitors.csv');
fprintf('--> Monitor File: %s\n',mon.fnm);
% _type of simulation_
mon.typ  = 'S';
fprintf('--> Type of Simulation: %s\n',mon.typ);
% _monitor identity_
mon.na = bhr.ns;
for m_ = 1:mon.na
    for n_ = 1:fnn.monn
        mon.(fni.monn{n_})(m_) = monn.(fni.monn{n_})(selected_case(m_));
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
fprintf('--> Components: ');
cellfun(@(x) fprintf('%s ',x),mon.rc);
fprintf('\n---------------------------------------------------------------\n');
mon.cp = mon.cp(:);
mon.ci = mon.ci(:);
mon.hyb = 'butter';

%% *MAP METADATA*
mon.map.flg = 1;
mon.map.typ = {'pga','pgv','psa'};
mon.map.stk = 95;
mon.map.vTn.psa = [0.1,0.2,0.3,0.5,0.75,1,2];
for i_=1:numel(mon.map.vTn.psa)
    fprintf('%.2f s',mon.map.vTn.psa(i_))
end
mon.map.vTn.rsd = 0.75;
mon.map.fnm = '/workdir/gattif/syn2ann_maps/map_mrn';

% _SP96 metadata_
mtd.sp96.na = mon.na;
for m_ = 1:mtd.sp96.na
    for n_ = 1:fnn.mtdd.sp96
        mtd.sp96.(fni.mtdd.sp96{n_})(m_,:) = ...
            mtdd.sp96.(fni.mtdd.sp96{n_})(selected_case(m_),:);
    end
end

save(sprintf('syn2ann_input_maps_%u.mat',NJB),'wd','mon','bhr','mtd','ann','hybrid_type','MAXIT');
