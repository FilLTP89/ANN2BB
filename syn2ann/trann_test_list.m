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
bhrr.st{1}.id = {'AMT'};
bhrr.st{1}.ni = {'';'AVD'};
bhrr.st{2}.id = {'NRC'};
bhrr.st{2}.ni = {'';'AVD'};
% _recorded events_
bhrr.st{1}.ev = {''};
bhrr.st{2}.ev = {''};
% _device list_
bhrr.st{1}.dv = {''};
bhrr.st{2}.dv = {''};
% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

%% *ANN METADATA ann*
ann.trn.nr = 6;
% gh <---> AB
ann.trn.mtd(1).TnC = 0.75;
ann.trn.mtd(1).cp = 'gh';
ann.trn.mtd(1).scl = 'AB';
% ud <---> AB
ann.trn.mtd(2).TnC = 0.75;
ann.trn.mtd(2).cp = 'ud';
ann.trn.mtd(2).scl = 'AB';
% gh <---> CD
ann.trn.mtd(3).TnC = 0.75;
ann.trn.mtd(3).cp = 'gh';
ann.trn.mtd(3).scl = 'CD';
% ud <---> CD
ann.trn.mtd(4).TnC = 0.75;
ann.trn.mtd(4).cp = 'ud';
ann.trn.mtd(4).scl = 'CD';
% gh <---> ALL
ann.trn.mtd(5).TnC = 0.75;
ann.trn.mtd(5).cp = 'gh';
ann.trn.mtd(5).scl = 'ALL';
% ud <---> ALL
ann.trn.mtd(6).TnC = 0.75;
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
tst.mtd.scl = {'ALL';'ALL';'ALL'};
% _corner period for each ANN_
tst.mtd.TnC = {0.75;0.75;0.75};
% _ANN motion component : gh,ud (geometric mean horizontal, updip)_
tst.mtd.cpn = {'gh';'gh';'ud'};
for i_ = 1:tst.mtd.nr
    tst.mtd.nl(i_) = {sprintf('net_%u_%s_%s.mat',...
        round(tst.mtd.TnC{i_}*100),tst.mtd.scl{i_},tst.mtd.cpn{i_})};
end