%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_case_list_: function to load all the metadata required to the 
% analyses to be run.
%% *N.B.*
% Need for:

%% *RECORDING STATION: bhrr*
% _station identity_

for i_=1:100
    bhrr.st{i_}.id = {'';''};
    bhrr.st{i_}.ni = {'';''};
    bhrr.st{i_}.ev = {'20120529.070002'};
    bhrr.st{i_}.dv = {''};
    bhrr.st{i_}.dv = {'itaca'};
end

% _bhrr field names_
fni.bhrr = fieldnames(bhrr);
fnn.bhrr = numel(fni.bhrr);

%% *MONITOR STATION: monn*
% _station identity_
monn.id = [16928,15045,1,2,18446,18437];

% _monn field names_
fni.monn = fieldnames(monn);
fnn.monn = numel(fni.monn);