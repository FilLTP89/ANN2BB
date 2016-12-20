%% *SIMULATION ON TRAINED ANN (NO PGV)*
fprintf('---------------------\n5. ANN\n---------------------\n');

% parse ann networks on each motion component
fprintf('--> Parsing\n');
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
