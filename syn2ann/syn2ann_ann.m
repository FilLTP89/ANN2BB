%% *TRAINED ANN*
fprintf('---------------------\n5. ANN\n---------------------\n');
% _ parse trained ANN network - horizontal directions_
% _parse trained ANN network - x direction_
% check motion components and load ann file
% list of the names of ann networks along each motion direction (according
% to mon.rs)

% parse ann networks on each motion component
fprintf('--> Parsing\n');
for j_ = 1:mon.nc
    [~,Locb]=ismember(mon.cp{j_},mon.rs);
    fn = fullfile(wd,'training',ann.mtd.nl{Locb});
    TnC  = ann.mtd.tc{Locb};
    ann.(mon.cp{j_}) = nn_parser(TnC,fn);
end

ann.cp = fieldnames(ann);
[~,~,ib] = intersect(hbs.mon.cp,ann.cp,'stable');
ann.cp = ann.cp(ib);

% _ apply trained ANN on hybrid accelerograms_
fprintf('--> Apply\n');
trs = ann2hbs_train_noPGV(hbs,ann);