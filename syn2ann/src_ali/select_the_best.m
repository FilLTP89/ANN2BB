%%% selects the best among 500.

sel_hyb_acc(:,1)=hyb_EW_acc(:,num_select1);
sel_hyb_acc(:,2)=hyb_NS_acc(:,num_select2);
sel_hyb_acc(:,3)=hyb_Z_acc(:,num_select3);

sel_hyb_PSA(:,1)=hyb_EW_PSA(:,num_select1);
sel_hyb_PSA(:,2)=hyb_NS_PSA(:,num_select2);
sel_hyb_PSA(:,3)=hyb_Z_PSA(:,num_select3);

% note that there are two different options, the above 6 lines are default,
% however, to represent more variability the below scheme may be used as
% well. Above, selects the best ones according to their scores in
% corresponding orientations (E,W,Z). Below, selects them according to
% overall combined performance.

% sel_hyb_acc(:,1)=hyb_EW_acc(:,num_select);
% sel_hyb_acc(:,2)=hyb_NS_acc(:,num_select);
% sel_hyb_acc(:,3)=hyb_Z_acc(:,num_select);
% 
% sel_hyb_PSA(:,1)=hyb_EW_PSA(:,num_select);
% sel_hyb_PSA(:,2)=hyb_NS_PSA(:,num_select);
% sel_hyb_PSA(:,3)=hyb_Z_PSA(:,num_select);

sel_hyb_vel(:,1)=cumsum(sel_hyb_acc(:,1)).*dt;
sel_hyb_vel(:,2)=cumsum(sel_hyb_acc(:,2)).*dt;
sel_hyb_vel(:,3)=cumsum(sel_hyb_acc(:,3)).*dt;

sel_hyb_dis(:,1)=cumsum(sel_hyb_vel(:,1)).*dt;
sel_hyb_dis(:,2)=cumsum(sel_hyb_vel(:,2)).*dt;
sel_hyb_dis(:,3)=cumsum(sel_hyb_vel(:,3)).*dt;

