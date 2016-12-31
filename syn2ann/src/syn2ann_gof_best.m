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
wgt = 1.00/mon.nc;
%% *LOOP OVER NUMBER OF MONITORS*
for i_=1:mon.na
    
    % FINAL CUMULATIVE GOF (ALL)
    err.gof{i_}.psa.all = zeros(MAXIT,1);
    
    %% *LOOP OVER NUMBER OF DIRECTIONS*
    for j_ = 1:mon.nc
        
        cpn = mon.cp{j_};
        
        %% *NORMALIZED AVERAGE ERROR OVER ITERATIONS*
        err.nav{i_}.psa.(cpn)(:,1) = err.avg{i_}.psa.(cpn)(:,1)./max(err.avg{i_}.psa.(cpn)(:,1));
        
        %% *NORMALIZED MAXIMUM ERROR OVER ITERATIONS*
        err.nmx{i_}.psa.(cpn)(:,1) = err.max{i_}.psa.(cpn)(:,1)./max(err.max{i_}.psa.(cpn)(:,1));
        
        %% *GOF OF SINGLE DIRECTION*
        %         err.gof{i_}.psa.(cpn)(:,1) = (2/3).*err.nav{i_}.psa.(cpn)(:,1)+...
        %             (1/3).*err.nmx{i_}.psa.(cpn)(:,1);
        err.gof{i_}.psa.(cpn)(:,1) = (1.00).*err.nav{i_}.psa.(cpn)(:,1)+...
            (0.00).*err.nmx{i_}.psa.(cpn)(:,1);
        
        %% *MINIMUM GOF OF SINGLE DIRECTION*
        [err.mgo{i_}.psa.(cpn),err.igo{i_}.psa.(cpn)] = min(err.gof{i_}.psa.(cpn)(:,1));
        
        
        %% *CUMULATIVE GOF OVER DIRECTIONS*
        % wgt = 0.25+0.25*(strcmpi(cpn,'ud'));
        err.gof{i_}.psa.all(:,1) = err.gof{i_}.psa.all(:,1) + ...
            wgt.*err.gof{i_}.psa.(cpn)(:,1);
        
    end
    
    %% *MINIMUM GOF OVER DIRECTIONS*
    [err.mgo{i_}.psa.all,err.igo{i_}.psa.all] = min(err.gof{i_}.psa.all(:,1));
    
end

hbs.bst = hbs.sps{err.igo{i_}.psa.all};
pbs.bst = pbs.sps{err.igo{i_}.psa.all};