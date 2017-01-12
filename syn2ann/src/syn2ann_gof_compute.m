%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_gof_compute_: function to compute gof structures.
%% *N.B.*
% Need for:

fprintf('============================\n');
fprintf('-------7. COMPUTE GOF-------\n');
fprintf('============================\n');

for i_=1:mon.na
    for j_ = 1:mon.nc
        cpn = mon.cp{j_};
        err.sng{i_}.psa.(cpn)(:,NIT) = ...
            abs(hbs.sps{NIT}.syn{i_}.psa.(cpn)([trs.sps.(cpn).tid;trs.sps.(cpn).iid])-...
            trs.sps.(cpn).syn{i_}.psa.(cpn)(:))./...
            trs.sps.(cpn).syn{i_}.psa.(cpn)(:);
        
        err.avg{i_}.psa.(cpn)(NIT,1) = mean(err.sng{i_}.psa.(cpn)(:,NIT));
        err.max{i_}.psa.(cpn)(NIT,1) = max(err.sng{i_}.psa.(cpn)(:,NIT));
    end
end