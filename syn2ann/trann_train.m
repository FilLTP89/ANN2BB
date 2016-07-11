%% *SET-UP*
fprintf('---------------------\n0. TRAIN ANN\n---------------------\n');
for i_ = 1:ann.nr
    train_ann_noPGV(ann.mtd(i_),ann.wd);
end