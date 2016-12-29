%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_gof_best_: function to compute the best gof.
%% *N.B.*
% Need for:

fprintf('============================\n');
fprintf('---------8. BEST GOF--------\n');
fprintf('============================\n');
for i_=1:mon.na
    % gof on all the directions
    err.gof{i_}.psa.all = zeros(MAXIT,1);
    for j_ = 1:mon.nc
        cpn = mon.cp{j_};
        
        % normalized average error over the iterations
        err.nav{i_}.psa.(cpn)(:,1) = err.avg{i_}.psa.(cpn)(:,1)./max(err.avg{i_}.psa.(cpn)(:,1));
        % normalized maximum error over the iterations
        err.nmx{i_}.psa.(cpn)(:,1) = err.max{i_}.psa.(cpn)(:,1)./max(err.max{i_}.psa.(cpn)(:,1));
        % gof along each direction
        err.gof{i_}.psa.(cpn)(:,1) = (2/3).*err.nav{i_}.psa.(cpn)(:,1)+...
            (1/3).*err.nmx{i_}.psa.(cpn)(:,1);
        % gof along each direction
        wgt = 0.25+0.25*(strcmpi(cpn,'ud'));
        
        err.gof{i_}.psa.all(:,1) = err.gof{i_}.psa.all(:,1) + ...
            wgt.*err.gof{i_}.psa.(cpn)(:,1);
        [err.mgo{i_}.psa.(cpn),err.igo{i_}.psa.(cpn)] = min(err.gof{i_}.psa.(cpn)(:,1));
    end
    [err.mgo{i_}.psa.all,err.igo{i_}.psa.all] = min(err.gof{i_}.psa.all(:,1));
    
end
