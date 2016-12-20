%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_ann_drive_: function to parse and simulate ANN upon numerical
% simulations
%% *N.B.*
% Need for:
% ismember.m,syn2ann_ann_parser.m,intersect.m,apply_ann2hbs_justPSA.m_
fprintf('---------------------\n5. SETUP ANN\n---------------------\n');

%% *PARSE ANN NETWORK ON EACH MOTION COMPONENT*
fprintf('--> Parsing \n');
for j_ = 1:mon.nc
    [~,ib] = ismember(mon.cp{j_},mon.rs);
    fnm = fullfile(wd,'training',ann.mtd.nl{ib});
    TnC  = ann.mtd.TnC{ib};
    ann.cp{j_} = mon.cp{j_};
    ann.(mon.cp{j_}) = syn2ann_ann_parser(TnC,fnm,mon.cp{j_},ann.mtd.scl{ib});
    ann.(mon.cp{j_}).nl = sprintf('%s_%u_%s',...
        ann.mtd.scl{ib},round(TnC*100),ann.cp{j_});
end
[~,~,ib] = intersect(pbs.org.mon.cp,ann.cp,'stable');
ann.cp   = ann.cp(ib);

%% *DEFINE TARGET SPECTRUM UPON NUMERICAL SIMULATIONS*
fprintf('--> Apply:\n');
trs.sps = apply_ann2hbs_justPSA(pbs.org,ann);