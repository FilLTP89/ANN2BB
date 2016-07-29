%% *SET-UP*
fprintf('---------------------\n0. TRAIN ANN (NO PGV-PGD)\n---------------------\n');
for i_ = 1:ann.nr
    train_ann_justPSA(ann.mtd(i_),ann.wd);
end