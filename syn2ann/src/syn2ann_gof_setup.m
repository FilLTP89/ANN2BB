%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_gof_setup_: function to setup gof structures.
%% *N.B.*
% Need for:

fprintf('============================\n');
fprintf('-------6. GOF SETUP-------\n');
fprintf('============================\n');

for i_=1:mon.na
    for j_ = 1:mon.nc
        cpn = mon.cp{j_};
        err.sng{i_}.tha.(cpn) = zeros(hbs.sps{1}.mon.ntm,MAXIT);
        err.sng{i_}.psa.(cpn) = zeros(hbs.sps{1}.mon.nT,MAXIT);
        err.avg{i_}.tha.(cpn) = zeros(MAXIT,1);
        err.avg{i_}.psa.(cpn) = zeros(MAXIT,1);
        err.max{i_}.tha.(cpn) = zeros(MAXIT,1);
        err.max{i_}.psa.(cpn) = zeros(MAXIT,1);
        err.nav{i_}.tha.(cpn) = zeros(MAXIT,1);
        err.nav{i_}.psa.(cpn) = zeros(MAXIT,1);
        err.nmx{i_}.tha.(cpn) = zeros(MAXIT,1);
        err.nmx{i_}.psa.(cpn) = zeros(MAXIT,1);
        gof{i_}.(cpn)         = zeros(MAXIT,1);
    end
    gof{i_}.all           = zeros(MAXIT,1);
end