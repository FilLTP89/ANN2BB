%% *SIMULATION ON TRAINED ANN (WITH PGV)*
fprintf('---------------------\n5. ANN (JUST PSA)\n---------------------\n');
% parse ann networks on each motion component
fprintf('--> Parsing (PGV included)\n');
ann_justPSA = ann;
for j_ = 1:mon.nc
    [~,Locb]=ismember(mon.cp{j_},mon.rs);
    fn = fullfile(wd,'training',ann_justPSA.mtd.nl.justPSA{Locb});
    TnC  = ann_justPSA.mtd.tc{Locb};
    ann_justPSA.(mon.cp{j_}) = nn_parser(TnC,fn);
end

ann_justPSA.cp = fieldnames(ann_justPSA);
[~,~,ib] = intersect(hbs.mon.cp,ann_justPSA.cp,'stable');
ann_justPSA.cp = ann_justPSA.cp(ib);

% _ apply trained ANN on hybrid accelerograms_
fprintf('--> Apply\n');
trs_justPSA = ann2hbs_train_justPSA(hbs,ann_justPSA);