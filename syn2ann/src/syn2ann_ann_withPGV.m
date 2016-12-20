%% *SIMULATION ON TRAINED ANN (WITH PGV)*
fprintf('---------------------\n5. ANN\n---------------------\n');

% parse ann networks on each motion component
fprintf('--> Parsing\n');
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
