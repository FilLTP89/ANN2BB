%% *LOAD TRAINED ANN*
fprintf('---------------------\n1. ANN (JUST PSA)\n---------------------\n');
% parse ann networks on each motion component
fprintf('--> Parsing \n');
for j_ = 1:tst.mtd.nr
    fn = fullfile(wd,'training',tst.mtd.nl{j_});
    TnC  = tst.mtd.TnC{j_};
    ann.tst{j_} = syn2ann_ann_parser(TnC,fn);
end