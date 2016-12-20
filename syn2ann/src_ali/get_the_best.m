%%% selects the best hybrid out of 50

if n_hyb==1
    hyb_EW_acc=zeros(length(hybrid.acc(:,1)),50);
    hyb_NS_acc=zeros(length(hybrid.acc(:,2)),50);
    hyb_Z_acc=zeros(length(hybrid.acc(:,3)),50);
    hyb_EW_PSA=zeros(length(hybrid.PSA(:,1)),50);
    hyb_NS_PSA=zeros(length(hybrid.PSA(:,2)),50);
    hyb_Z_PSA=zeros(length(hybrid.PSA(:,3)),50);
    err_EW_PSA=zeros(length(hybrid.PSA(:,1)),50);
    err_NS_PSA=zeros(length(hybrid.PSA(:,2)),50);
    err_Z_PSA=zeros(length(hybrid.PSA(:,3)),50);
    max_err_EW=zeros(1,50);mean_err_EW=zeros(1,50);
    max_err_NS=zeros(1,50);mean_err_NS=zeros(1,50);
    max_err_Z=zeros(1,50);mean_err_Z=zeros(1,50);
    nor_max_err_EW=zeros(1,50);nor_mean_err_EW=zeros(1,50);
    nor_max_err_NS=zeros(1,50);nor_mean_err_NS=zeros(1,50);
    nor_max_err_Z=zeros(1,50);nor_mean_err_Z=zeros(1,50);
    score_EW=zeros(1,50);score_NS=zeros(1,50);score_Z=zeros(1,50);
    score_ALL=zeros(1,50);
end


hyb_EW_acc(:,n_hyb)=hybrid.acc(:,1);
hyb_NS_acc(:,n_hyb)=hybrid.acc(:,2);
hyb_Z_acc(:,n_hyb)=hybrid.acc(:,3);

hyb_EW_PSA(:,n_hyb)=hybrid.PSA(:,1);
hyb_NS_PSA(:,n_hyb)=hybrid.PSA(:,2);
hyb_Z_PSA(:,n_hyb)=hybrid.PSA(:,3);

if n_hyb==50
  for k=1:50
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
        if k==50
            for m=1:50
            nor_max_err_EW(m)=max_err_EW(m)/max(max_err_EW(1:50));
            nor_max_err_NS(m)=max_err_NS(m)/max(max_err_NS(1:50));
            nor_max_err_Z(m)=max_err_Z(m)/max(max_err_Z(1:50));
            nor_mean_err_EW(m)=mean_err_EW(m)/max(mean_err_EW(1:50));
            nor_mean_err_NS(m)=mean_err_NS(m)/max(mean_err_NS(1:50));
            nor_mean_err_Z(m)=mean_err_Z(m)/max(mean_err_Z(1:50));
            score_EW(m)=0.33*nor_max_err_EW(m)+0.67*nor_mean_err_EW(m);
            score_NS(m)=0.33*nor_max_err_NS(m)+0.67*nor_mean_err_NS(m);
            score_Z(m)=0.33*nor_max_err_Z(m)+0.67*nor_mean_err_Z(m);
            score_ALL(m)=0.50*score_Z(m)+0.25*score_EW(m)+0.25*score_NS(m);
            end
            [min_select,num_select]=min(score_ALL);
        end
  end
end

sel_hyb_EW_acc=hyb_EW_acc(:,num_select);
sel_hyb_NS_acc=hyb_NS_acc(:,num_select);
sel_hyb_Z_acc=hyb_Z_acc(:,num_select);

sel_hyb_EW_PSA=hyb_EW_PSA(:,num_select);
sel_hyb_NS_PSA=hyb_NS_PSA(:,num_select);
sel_hyb_Z_PSA=hyb_Z_PSA(:,num_select);