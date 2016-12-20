
    hyb_EW_acc=zeros(length(hybrid.acc(:,1)),500);
    hyb_NS_acc=zeros(length(hybrid.acc(:,2)),500);
    hyb_Z_acc=zeros(length(hybrid.acc(:,3)),500);
    hyb_EW_PSA=zeros(length(hybrid.PSA(:,1)),500);
    hyb_NS_PSA=zeros(length(hybrid.PSA(:,2)),500);
    hyb_Z_PSA=zeros(length(hybrid.PSA(:,3)),500);
    err_EW_PSA=zeros(length(hybrid.PSA(:,1)),500);
    err_NS_PSA=zeros(length(hybrid.PSA(:,2)),500);
    err_Z_PSA=zeros(length(hybrid.PSA(:,3)),500);
    max_err_EW=zeros(1,50);mean_err_EW=zeros(1,500);
    max_err_NS=zeros(1,50);mean_err_NS=zeros(1,500);
    max_err_Z=zeros(1,50);mean_err_Z=zeros(1,500);
    nor_max_err_EW=zeros(1,50);nor_mean_err_EW=zeros(1,500);
    nor_max_err_NS=zeros(1,50);nor_mean_err_NS=zeros(1,500);
    nor_max_err_Z=zeros(1,50);nor_mean_err_Z=zeros(1,500);
    score_EW=zeros(1,50);score_NS=zeros(1,50);score_Z=zeros(1,500);
    score_ALL=zeros(1,500);
    
    T =[tar_Tn inp_Tn];
