%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_test_list_: function to load all the metadata required to the 
% analyses to be run.
%% *N.B.*
% Need for:

%% *RECORDING STATION: bhrr*
% _station identity_
bhrr.st{1}.id = {'ACC'};
bhrr.st{1}.ni = {'IT';'HG'};
bhrr.st{2}.id = {'ACC'};
bhrr.st{2}.ni = {'IT';'HG'};
bhrr.st{3}.id = {'ACC'};
bhrr.st{3}.ni = {'IT';'HG'};
bhrr.st{4}.id = {'NRC'};
bhrr.st{4}.ni = {'IT';'HG'};
bhrr.st{5}.id = {'NRC'};
bhrr.st{5}.ni = {'IT';'HG'};
bhrr.st{6}.id = {'NRC'};
bhrr.st{6}.ni = {'IT';'HG'};
bhrr.st{7}.id = {'NRC'};
bhrr.st{7}.ni = {'IT';'HG'};
bhrr.st{8}.id = {'NRC'};
bhrr.st{8}.ni = {'IT';'HG'};
bhrr.st{9}.id = {'KMMH16'};
bhrr.st{9}.ni = {'';''};
bhrr.st{10}.id = {'KMMH16'};
bhrr.st{10}.ni = {'';''};
bhrr.st{11}.id = {'KMMH16'};
bhrr.st{11}.ni = {'';''};
bhrr.st{12}.id = {'KMMH16'};
bhrr.st{12}.ni = {'';''};
bhrr.st{13}.id = {'KMMH16'};
bhrr.st{13}.ni = {'';''};
bhrr.st{14}.id = {'KMMH14'};
bhrr.st{14}.ni = {'';''};
bhrr.st{15}.id = {'KMMH14'};
bhrr.st{15}.ni = {'';''};
bhrr.st{16}.id = {'KMMH14'};
bhrr.st{16}.ni = {'';''};
bhrr.st{17}.id = {'KMMH16'};
bhrr.st{17}.ni = {'';''};
bhrr.st{18}.id = {'KMMH16'};
bhrr.st{18}.ni = {'';''};
bhrr.st{19}.id = {'KMMH16'};
bhrr.st{19}.ni = {'';''};
bhrr.st{20}.id = {'KMMH16'};
bhrr.st{20}.ni = {'';''};
bhrr.st{21}.id = {'KMMH16'};
bhrr.st{21}.ni = {'';''};
bhrr.st{22}.id = {'KMMH14'};
bhrr.st{22}.ni = {'';''};
bhrr.st{23}.id = {'KMMH14'};
bhrr.st{23}.ni = {'';''};
bhrr.st{24}.id = {'KMMH14'};
bhrr.st{24}.ni = {'';''};
% _recorded events_
bhrr.st{1}.ev  = {'20161026.171036'};
bhrr.st{2}.ev  = {'20161026.191806'};
bhrr.st{3}.ev  = {'20161030.064018'};
bhrr.st{4}.ev  = {'20160824.013632'};
bhrr.st{5}.ev  = {'20160824.023329'};
bhrr.st{6}.ev  = {'20161026.171036'};
bhrr.st{7}.ev  = {'20161026.191806'};
bhrr.st{8}.ev  = {'20161030.064018'};
bhrr.st{9}.ev  = {'1604142126'};
bhrr.st{10}.ev = {'1604142207'};
bhrr.st{11}.ev = {'1604150003'};
bhrr.st{12}.ev = {'1604160125'};
bhrr.st{13}.ev = {'1604160146'};
bhrr.st{14}.ev = {'1604142126'};
bhrr.st{15}.ev = {'1604150003'};
bhrr.st{16}.ev = {'1604160125'};
bhrr.st{17}.ev = {'1604142126'};
bhrr.st{18}.ev = {'1604142207'};
bhrr.st{19}.ev = {'1604150003'};
bhrr.st{20}.ev = {'1604160125'};
bhrr.st{21}.ev = {'1604160146'};
bhrr.st{22}.ev = {'1604142126'};
bhrr.st{23}.ev = {'1604150003'};
bhrr.st{24}.ev = {'1604160125'};
% _device list_
bhrr.st{1}.dv = {''};
bhrr.st{2}.dv = {''};
bhrr.st{3}.dv = {''};
bhrr.st{4}.dv = {''};
bhrr.st{5}.dv = {''};
bhrr.st{6}.dv = {''};
bhrr.st{7}.dv = {''};
bhrr.st{8}.dv = {''};
bhrr.st{9}.dv  = {'1'};
bhrr.st{10}.dv = {'1'};
bhrr.st{11}.dv = {'1'};
bhrr.st{12}.dv = {'1'};
bhrr.st{13}.dv = {'1'};
bhrr.st{14}.dv = {'1'};
bhrr.st{15}.dv = {'1'};
bhrr.st{16}.dv = {'1'};
bhrr.st{17}.dv = {'2'};
bhrr.st{18}.dv = {'2'};
bhrr.st{19}.dv = {'2'};
bhrr.st{20}.dv = {'2'};
bhrr.st{21}.dv = {'2'};
bhrr.st{22}.dv = {'2'};
bhrr.st{23}.dv = {'2'};
bhrr.st{24}.dv = {'2'};
% _database_
bhrr.st{1}.tp = {'itaca'};
bhrr.st{2}.tp = {'itaca'};
bhrr.st{3}.tp = {'itaca'};
bhrr.st{4}.tp = {'itaca'};
bhrr.st{5}.tp = {'itaca'};
bhrr.st{6}.tp = {'itaca'};
bhrr.st{7}.tp = {'itaca'};
bhrr.st{8}.tp = {'itaca'};
bhrr.st{9}.tp = {'kiknet'};
bhrr.st{10}.tp = {'kiknet'};
bhrr.st{11}.tp = {'kiknet'};
bhrr.st{12}.tp = {'kiknet'};
bhrr.st{13}.tp = {'kiknet'};
bhrr.st{14}.tp = {'kiknet'};
bhrr.st{15}.tp = {'kiknet'};
bhrr.st{16}.tp = {'kiknet'};
bhrr.st{17}.tp = {'kiknet'};
bhrr.st{18}.tp = {'kiknet'};
bhrr.st{19}.tp = {'kiknet'};
bhrr.st{20}.tp = {'kiknet'};
bhrr.st{21}.tp = {'kiknet'};
bhrr.st{22}.tp = {'kiknet'};
bhrr.st{23}.tp = {'kiknet'};
bhrr.st{24}.tp = {'kiknet'};
% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

%% *ANN METADATA ann*
ann.trn.nr = 6;
% gh <---> AB
ann.trn.mtd(1).TnC = 0.5;
ann.trn.mtd(1).cp = 'gh';
ann.trn.mtd(1).scl = 'AB';
% ud <---> AB
ann.trn.mtd(2).TnC = 0.5;
ann.trn.mtd(2).cp = 'ud';
ann.trn.mtd(2).scl = 'AB';
% gh <---> CD
ann.trn.mtd(3).TnC = 0.5;
ann.trn.mtd(3).cp = 'gh';
ann.trn.mtd(3).scl = 'CD';
% ud <---> CD
ann.trn.mtd(4).TnC = 0.5;
ann.trn.mtd(4).cp = 'ud';
ann.trn.mtd(4).scl = 'CD';
% gh <---> ALL
ann.trn.mtd(5).TnC = 0.5;
ann.trn.mtd(5).cp = 'gh';
ann.trn.mtd(5).scl = 'ALL';
% ud <---> ALL
ann.trn.mtd(6).TnC = 0.5;
ann.trn.mtd(6).cp = 'ud';
ann.trn.mtd(6).scl = 'ALL';
% _database_
for i_ = 1:ann.trn.nr
    ann.trn.mtd(i_).dbn = dbn;
end

%% *DEFINE ANN TEST METADATA*
% _number of ann to be tested_
tst.mtd.nr = 3;
% _site class considered : ALL,AB,CD_
tst.mtd.scl = {'AB';'AB';'AB'};
% _corner period for each ANN_
tst.mtd.TnC = {0.5;0.75;1};
% _ANN motion component : gh,ud (geometric mean horizontal, updip)_
tst.mtd.cpp = {'ud';'ud';'ud'};
for i_ = 1:tst.mtd.nr
    tst.mtd.nl(i_) = {sprintf('net_%u_%s_%s_new.mat',...
        round(tst.mtd.TnC{i_}*100),tst.mtd.scl{i_},tst.mtd.cpp{i_})};
end
