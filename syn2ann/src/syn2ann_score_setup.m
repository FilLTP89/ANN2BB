%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_score_setup_: function to setup gof structures.
%% *N.B.*
% Need for:

fprintf('============================\n');
fprintf('-------6. SCORE SETUP-------\n');
fprintf('============================\n');

for i_=1:mon.na
    for j_ = 1:mon.nc
        cpn = mon.cp{j_};
        err.sng{i_}.tha.(cpn) = zeros(hbs.sps.mon.ntm,MAXIT);
        err.sng{i_}.psa.(cpn) = zeros(hbs.sps.mon.nT,MAXIT);
        err.avg{i_}.tha.(cpn) = zeros(1,MAXIT);
        err.avg{i_}.psa.(cpn) = zeros(1,MAXIT);
        err.max{i_}.tha.(cpn) = zeros(1,MAXIT);
        err.max{i_}.psa.(cpn) = zeros(1,MAXIT);
        err.nav{i_}.tha.(cpn) = zeros(1,MAXIT);
        err.nav{i_}.psa.(cpn) = zeros(1,MAXIT);
        err.nmx{i_}.tha.(cpn) = zeros(1,MAXIT);
        err.nmx{i_}.psa.(cpn) = zeros(1,MAXIT);
        gof{i_}.(cpn)         = zeros(1,MAXIT);
    end
    gof{i_}.all               = zeros(1,MAXIT);
end