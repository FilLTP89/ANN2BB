%%% scores the 500 hybrid motion according to their agreement with target
%%% spectrum
n_hyb=x;

hyb_EW_acc(:,n_hyb)=hybrid.acc(:,1);
hyb_NS_acc(:,n_hyb)=hybrid.acc(:,2);
hyb_Z_acc(:,n_hyb)=hybrid.acc(:,3);

hyb_EW_PSA(:,n_hyb)=hybrid.PSA(:,1)*100; %cm/s2
hyb_NS_PSA(:,n_hyb)=hybrid.PSA(:,2)*100;
hyb_Z_PSA(:,n_hyb)=hybrid.PSA(:,3)*100;

if n_hyb==500
  for k=1:500
      for j=1:length(T)
      err_EW_PSA(j,k)=abs(hyb_EW_PSA(j,k)-PSA_target(j,1))/PSA_target(j,1);
      err_NS_PSA(j,k)=abs(hyb_NS_PSA(j,k)-PSA_target(j,2))/PSA_target(j,2);
      err_Z_PSA(j,k)=abs(hyb_Z_PSA(j,k)-PSA_target(j,3))/PSA_target(j,3);
      end
      max_err_EW(k)=max(err_EW_PSA(:,k));
      max_err_NS(k)=max(err_NS_PSA(:,k));
      max_err_Z(k)=max(err_Z_PSA(:,k));
      mean_err_EW(k)=mean(err_EW_PSA(:,k));
      mean_err_NS(k)=mean(err_NS_PSA(:,k));
      mean_err_Z(k)=mean(err_Z_PSA(:,k));
        if k==500
            for m=1:500
            nor_max_err_EW(m)=max_err_EW(m)/max(max_err_EW(1:500));
            nor_max_err_NS(m)=max_err_NS(m)/max(max_err_NS(1:500));
            nor_max_err_Z(m)=max_err_Z(m)/max(max_err_Z(1:500));
            nor_mean_err_EW(m)=mean_err_EW(m)/max(mean_err_EW(1:500));
            nor_mean_err_NS(m)=mean_err_NS(m)/max(mean_err_NS(1:500));
            nor_mean_err_Z(m)=mean_err_Z(m)/max(mean_err_Z(1:500));
            score_EW(m)=0.0*nor_max_err_EW(m)+1.00*nor_mean_err_EW(m);
            score_NS(m)=0.0*nor_max_err_NS(m)+1.00*nor_mean_err_NS(m);
            score_Z(m)=0.0*nor_max_err_Z(m)+1.00*nor_mean_err_Z(m);
            score_ALL(m)=0.33*score_Z(m)+0.33*score_EW(m)+0.33*score_NS(m);
            end
%             [min_select,num_select]=min(score_ALL);
            [min_select,num_select1]=min(score_EW);
            [min_select,num_select2]=min(score_NS);
            [min_select,num_select3]=min(score_Z);
            [min_select,num_select]=min(score_ALL);
        end
  end
end


