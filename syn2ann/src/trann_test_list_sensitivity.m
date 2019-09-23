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
    
stss = {'G07','G08','G09','G10','1G1','SG1','SG2','SG3','SG4','NIG016','NIG018','5G1'};

for i_=1:200
    bhrr.st{i_}.id = stss(5);
    bhrr.st{i_}.ni = {'';''};

    % _recorded events_
    bhrr.st{i_}.ev  = {''};
    % _device list_
    bhrr.st{i_}.dv = {''};
    % _database_
    bhrr.st{i_}.tp = {''};
    % psa variation
    %% [MODIFY] FLUCTUATION
    bhrr.sns{i_}.vTn=mat.T_ann;
    bhrr.sns{i_}.mus=mat.mus(:,1);
    bhrr.sns{i_}.psa=mat.psa_new(:,i_);
    bhrr.sns{i_}.sig=mat.sigmas;
end

monn.id =[(28:31)';28;(19:22)';14;16;23];
% _monn field names_
fni.monn = fieldnames(monn);
fnn.monn = numel(fni.monn);
% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

%% *DEFINE ANN TEST METADATA (CUSTOMIZE)*
% _number of ann to be tested_
tst.mtd.nr = 2;
% _type of testing (scl: by site class;TnC: by natural period)_
% _NB: when testing by TnC, be careful to be coherent with ANN (same 
% characteristics, but different TnC)_
tst.typ_cmp = 'TnC';
% select ann metadata to test
if strcmpi(tst.typ_cmp,'TnC')
    % _select site class considered : ALL,AB,CD_
    tst.mtd.scl = {'ALL';'ALL'};
    % _select corner period for each ANN_
    tst.mtd.TnC = {0.75;0.25};
    % _select ANN motion component : gh,ud (geometric mean horizontal, updip)_
    tst.mtd.cpp = {'gh';'gh'};
    tst.mtd.Tno = [0.75;0.25];
elseif strcmpi(tst.typ_cmp,'scl')
    tst.mtd.scl = {'ALL';'AB';'CD'};
    % _corner period for each ANN_
    tst.mtd.TnC = {0.75;0.75;0.75};
    % _ANN motion component : gh,ud (geometric mean horizontal, updip)_
    tst.mtd.cpp = {'ud';'ud';'ud'};
    tst.mtd.Tno = [0.75;0.75;0.75];
end
% for i_ = 1:tst.mtd.nr
%     tst.mtd.nl(i_) = {sprintf('net_%u_%s_%s_new.mat',...
%         round(tst.mtd.TnC{i_}*100),tst.mtd.scl{i_},tst.mtd.cpp{i_})};
% end
for i_ = 1:tst.mtd.nr
    tst.mtd.nl(i_) = {sprintf('net_%u_%s_%s_30n.mat',...
    round(tst.mtd.TnC{i_}*100),tst.mtd.scl{i_},tst.mtd.cpp{i_})};
end
