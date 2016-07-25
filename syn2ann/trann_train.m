%% *SET-UP*
fprintf('---------------------\n0. TRAIN ANN (WITH PGV)\n---------------------\n');
for i_ = 1:ann.nr
    train_ann_withPGV(ann.mtd(i_),ann.wd);
end

fprintf('---------------------\n0. TRAIN ANN (NO PGV)\n---------------------\n');
for i_ = 1:ann.nr
    train_ann_noPGV(ann.mtd(i_),ann.wd);
end