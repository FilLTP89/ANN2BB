%% *SIMULATION ON TRAINED ANN (WITH PGV)*
fprintf('---------------------\n5. ANN\n---------------------\n');

% parse ann networks on each motion component
fprintf('--> Parsing (PGV included)\n');
ann_withPGV = ann;
for j_ = 1:mon.nc
    [~,Locb]=ismember(mon.cp{j_},mon.rs);
    fn = fullfile(wd,'training',ann_withPGV.mtd.nl.withPGV{Locb});
    TnC  = ann_withPGV.mtd.tc{Locb};
    ann_withPGV.(mon.cp{j_}) = nn_parser(TnC,fn);
end

ann_withPGV.cp = fieldnames(ann_withPGV);
[~,~,ib] = intersect(hbs.mon.cp,ann_withPGV.cp,'stable');
ann_withPGV.cp = ann_withPGV.cp(ib);

% _ apply trained ANN on hybrid accelerograms_
fprintf('--> Apply\n');
trs_withPGV = ann2hbs_train_withPGV(hbs,ann_withPGV);

%% *SIMULATION ON TRAINED ANN (NO PGV)*
% parse ann networks on each motion component
fprintf('--> Parsing (no PGV included)\n');
ann_noPGV = ann;
for j_ = 1:mon.nc
    [~,Locb]=ismember(mon.cp{j_},mon.rs);
    fn = fullfile(wd,'training',ann_noPGV.mtd.nl.noPGV{Locb});
    TnC  = ann_noPGV.mtd.tc{Locb};
    ann_noPGV.(mon.cp{j_}) = nn_parser(TnC,fn);
end

ann_noPGV.cp = fieldnames(ann_noPGV);
[~,~,ib] = intersect(hbs.mon.cp,ann_noPGV.cp,'stable');
ann_noPGV.cp = ann_noPGV.cp(ib);

% _ apply trained ANN on hybrid accelerograms_
fprintf('--> Apply\n');
trs_noPGV = ann2hbs_train_noPGV(hbs,ann_noPGV);

