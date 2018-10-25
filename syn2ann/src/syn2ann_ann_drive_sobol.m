%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_ann_drive_: function to parse and simulate ANN upon PBS numerical
% simulations
%% *N.B.*
% Need for:
% ismember.m,syn2ann_ann_parser.m,intersect.m,apply_ann2hbs_justPSA.m_
fprintf('============================\n');
fprintf('------3. ANN SIMULATIONS----\n');
fprintf('============================\n');

%% *PARSE ANN NETWORK ON EACH MOTION COMPONENT*
fprintf('--> Parsing ANN \n');
for j_ = 1:mon.nc
    % find direction according to reference system mon.rs
    [~,ib] = ismember(mon.cp{j_},mon.rs);
    fnm = fullfile(wd,'training',ann.mtd.nl{ib});
    
    TnC  = ann.mtd.TnC{ib};
    ann.cp{j_} = mon.cp{j_};
    ann.(mon.cp{j_}) = syn2ann_ann_parser(TnC,fnm,mon.cp{j_},ann.mtd.scl{ib});
    ann.(mon.cp{j_}).nl = sprintf('%s_%u_%s',...
        ann.mtd.scl{ib},round(TnC*100),ann.cp{j_});
    ann.(mon.cp{j_}).tol = ann.mtd.tol(ib);
    ann.(mon.cp{j_}).nit = ann.mtd.nit(ib);
end
% [~,~,ib] = intersect(pbs.org.mon.cp,ann.cp,'stable');
% ann.cp   = ann.cp(ib);

%% *DEFINE TARGET SPECTRUM UPON NUMERICAL SIMULATIONS*
fprintf('--> Apply ANN:\n');

pdf = 1;
cv = 10.0/100.0;
[trs.sps,VVarEntree] = prepare_ann2hbs_sobol(pbs.org,ann,pdf,cv);
