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
    [~,err.mon{i_}.idx,~] = intersect(hbs.sps.mon.vTn{i_},trs.sps.(cpn).mon.vTn);
    for j_ = 1:mon.nc
        cpn = mon.cp{j_};
        
        err.sng{i_}.psa.(cpn)(:,NIT) = ...
            abs(hbs.sps.syn{i_}.psa.(cpn)(:)-...
            trs.sps.(cpn).syn{i_}.psa.(cpn)(:))./...
            trs.sps.(cpn).syn{i_}.psa.(cpn)(:);
        
        err.avg{i_}.psa.(cpn)(NIT,1) = mean(err.sng{i_}.psa.(cpn)(:,NIT));
        err.max{i_}.psa.(cpn)(NIT,1) = max(err.sng{i_}.psa.(cpn)(:,NIT));
        err.nav{i_}.psa.(cpn)(NIT,1) = zeros(MAXIT,1);
        err.nmx{i_}.psa.(cpn)(NIT,1) = zeros(MAXIT,1);
        gof{i_}.(cpn)         = zeros(MAXIT,1);
    end
    gof{i_}.all               = zeros(MAXIT,1);
end

% for j=1:length(T)
%     err_EW_PSA(j,k)=abs(hyb_EW_PSA(j,k)-PSA_target(j,1))/PSA_target(j,1);
%     err_NS_PSA(j,k)=abs(hyb_NS_PSA(j,k)-PSA_target(j,2))/PSA_target(j,2);
%     err_Z_PSA(j,k)=abs(hyb_Z_PSA(j,k)-PSA_target(j,3))/PSA_target(j,3);
% end
% max_err_EW(k)=max(err_EW_PSA(:,k));
% max_err_NS(k)=max(err_NS_PSA(:,k));
% max_err_Z(k)=max(err_Z_PSA(:,k));
% mean_err_EW(k)=mean(err_EW_PSA(:,k));
% mean_err_NS(k)=mean(err_NS_PSA(:,k));
% mean_err_Z(k)=mean(err_Z_PSA(:,k));
% if k==MAXIT
%     for m=1:MAXIT
%         nor_max_err_EW(m)=max_err_EW(m)/max(max_err_EW(1:MAXIT));
%         nor_max_err_NS(m)=max_err_NS(m)/max(max_err_NS(1:MAXIT));
%         nor_max_err_Z(m)=max_err_Z(m)/max(max_err_Z(1:MAXIT));
%         nor_mean_err_EW(m)=mean_err_EW(m)/max(mean_err_EW(1:MAXIT));
%         nor_mean_err_NS(m)=mean_err_NS(m)/max(mean_err_NS(1:MAXIT));
%         nor_mean_err_Z(m)=mean_err_Z(m)/max(mean_err_Z(1:MAXIT));
%         score_EW(m)=0.0*nor_max_err_EW(m)+1.00*nor_mean_err_EW(m);
%         score_NS(m)=0.0*nor_max_err_NS(m)+1.00*nor_mean_err_NS(m);
%         score_Z(m)=0.0*nor_max_err_Z(m)+1.00*nor_mean_err_Z(m);
%         score_ALL(m)=0.33*score_Z(m)+0.33*score_EW(m)+0.33*score_NS(m);
%     end
%     %             [min_select,num_select]=min(score_ALL);
%     [min_select,num_select1]=min(score_EW);
%     [min_select,num_select2]=min(score_NS);
%     [min_select,num_select3]=min(score_Z);
%     [min_select,num_select]=min(score_ALL);
% end
% end
% end