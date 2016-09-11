NS = 1;

% _station identity_
bhrr.st{1}.id = {'MRN'};
bhrr.st{1}.ni = {'CIT';'HN'};
bhrr.st{2}.id = {'MIR'};
bhrr.st{2}.ni = {'TV';'HN'};
bhrr.st{3}.id = {'AQK'};
bhrr.st{3}.ni = {'IT';'HN'};
bhrr.st{4}.id = {'AQU'};
bhrr.st{4}.ni = {'MN';'HL'};

% _recorded events_
bhrr.st{1}.ev = {'20120529.070002'};
bhrr.st{2}.ev = {'20120529.070002'};
bhrr.st{3}.ev = {'20090406.013240'};
bhrr.st{4}.ev = {'20090406.013240'};

% _device list_
bhrr.st{1}.dv = {''};
bhrr.st{2}.dv = {'08'};
bhrr.st{3}.dv = {''};
bhrr.st{4}.dv = {''};


fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

monn.id = [16928,15045,1,2];
monn.fa = [1.3,1.3,1.3,1.3]; % in Hz
monn.fb = [1.7,1.5,1.5,1.5]; % in Hz

fni.monn = fieldnames(monn);
fnn.monn = numel(fni.monn);

mtdd.mw = [6,6,6.3,6.3];
mtdd.dtm_sp96 = [0.01,0.01,0.01,0.01];
mtdd.scc = [2,1,2,1];
mtdd.sst = zeros(4,1);
mtdd.scl = 0.01*ones(4,1);

fni.mtdd = fieldnames(mtdd);
fnn.mtdd = numel(fni.mtdd);