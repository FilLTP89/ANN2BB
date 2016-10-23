%% *SIMULATION ON TRAINED ANN (WITH PGV)*
fprintf('---------------------\n5. ANN (JUST PSA)\n---------------------\n');
% parse ann networks on each motion component
fprintf('--> Parsing \n');
for j_ = 1:mon.nc
    [~,ib] = ismember(mon.cp{j_},mon.rs);
    fn = fullfile(wd,'training',ann.mtd.nl{ib});
    TnC  = ann.mtd.TnC{ib};
    ann.(mon.cp{j_}) = syn2ann_ann_parser(TnC,fn);
    disp('CHECK DIRECTIONS')
    keyboard
end

ann.cp   = fieldnames(ann);

switch lower(hybrid_type)
    case 'sp96'
        %
        % _SABETTA & PUGLIESE 1996_
        %
        [~,~,ib] = intersect(hbs.sps.mon.cp,ann.cp,'stable');
        ann.cp   = ann.cp(ib);
        % _ apply trained ANN on hybrid accelerograms_
        fprintf('--> Apply\n');
        trs.sps = ann2hbs_train_justPSA(hbs.sps,ann);
    case 'exsim'
        %
        % _EXSIM_
        %
        [~,~,ib] = intersect(hbs.exs.mon.cp,ann.cp,'stable');
        ann.cp   = ann.cp(ib);
        % _ apply trained ANN on hybrid accelerograms_
        fprintf('--> Apply\n');
        trs.exs = ann2hbs_train_justPSA(hbs.exs,ann);
    case 'both'
        %
        % _SABETTA & PUGLIESE 1996_
        %
        [~,~,ib] = intersect(hbs.sps.mon.cp,ann.cp,'stable');
        ann.cp   = ann.cp(ib);
        % _ apply trained ANN on hybrid accelerograms_
        fprintf('--> Apply\n');
        trs.sps = ann2hbs_train_justPSA(hbs.sps,ann);
        %
        % _EXSIM_
        %
        [~,~,ib] = intersect(hbs.exs.mon.cp,ann.cp,'stable');
        ann.cp   = ann.cp(ib);
        % _ apply trained ANN on hybrid accelerograms_
        fprintf('--> Apply\n');
        trs.exs = ann2hbs_train_justPSA(hbs.exs,ann);
end