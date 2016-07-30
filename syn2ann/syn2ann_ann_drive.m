%% *SIMULATION ON TRAINED ANN (WITH PGV)*
fprintf('---------------------\n5. ANN (JUST PSA)\n---------------------\n');
% parse ann networks on each motion component
fprintf('--> Parsing \n');
for j_ = 1:mon.nc
    [~,ib] = ismember(mon.cp{j_},mon.rs);
    fn = fullfile(wd,'training',ann.mtd.nl{ib});
    TnC  = ann.mtd.TnC{ib};
    ann.(mon.cp{j_}) = syn2ann_ann_parser(TnC,fn);
end

ann.cp   = fieldnames(ann);
[~,~,ib] = intersect(hbs.mon.cp,ann.cp,'stable');
ann.cp   = ann.cp(ib);

% _ apply trained ANN on hybrid accelerograms_
fprintf('--> Apply\n');
trs = ann2hbs_train_justPSA(hbs,ann);