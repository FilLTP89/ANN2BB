%% * COMPARISON: Numerical Simulation vs Real Record*

flag = 0;
flag_TH = 0;
flag_RS = 0;
if cfr_record == 1
    if flag
       i_fig = 1; 
       for i = 1:3 % e,n,z components 
           %% **Comparison: TIME HISTORIES**
           if flag_TH
               %% ***Comparison: ACCELERATION TH***
               ymax_acc = max(max(abs(num_sim.acc_orig(:,i))),max(abs(record.acc(:,i))));
               % Single plots
               figure(i_fig*i)
               subplot(2,1,1)
               plot(num_sim.t_vec,num_sim.acc_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('a(t) [m/s^2]');
               xlim([0,40]); 
               ylim([-ymax_acc,ymax_acc]);
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
               subplot(2,1,2)
               plot(record.t_vec,record.acc(:,i),'r','Linewidth',1);
               legendInfo{1} = ('record');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('a(t) [m/s^2]');
               xlim([0,40]);
               ylim([-ymax_acc,ymax_acc]);
               legend(legendInfo);
               % Comparison plot
               figure(i_fig*i+5)
               plot(num_sim.t_vec,num_sim.acc_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               hold on 
               plot(record.t_vec,record.acc(:,i),'r','Linewidth',1);
               legendInfo{2} = ('record');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('a(t) [m/s^2]');
               ylim([-ymax_acc,ymax_acc]);
               xlim([0,40]);
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
               %% ***Comparison: VELOCITY TH***
               ymax_vel = max(max(abs(num_sim.vel_orig(:,i))),max(abs(record.vel(:,i))));
               % Single plots
               figure(i_fig*i+10)
               subplot(2,1,1)
               plot(num_sim.t_vec,num_sim.vel_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('v(t) [m/s]');
               ylim([-ymax_vel,ymax_vel]);
               xlim([0,40]);
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
               subplot(2,1,2)
               plot(record.t_vec,record.vel(:,i),'r','Linewidth',1);
               legendInfo{1} = ('record');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('v(t) [m/s]');
               xlim([0,40]);
               ylim([-ymax_vel,ymax_vel]);
               legend(legendInfo);
               % Comparison plot
               figure(i_fig*i+10+5)
               plot(num_sim.t_vec,num_sim.vel_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               hold on 
               plot(record.t_vec,record.vel(:,i),'r','Linewidth',1);
               legendInfo{2} = ('record');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('v(t) [m/s]');
               ylim([-ymax_vel,ymax_vel]);
               xlim([0,40]);
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
               %% ***Comparison: DISPLACEMNT TH***
               ymax_dis = max(max(abs(num_sim.dis_orig(:,i))),max(abs(record.dis(:,i))));
               % Single plots
               figure(i_fig*i+20)
               subplot(2,1,1)
               plot(num_sim.t_vec,num_sim.dis_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('d(t) [m]');
               xlim([0,40]);
               ylim([-ymax_dis,ymax_dis])
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
               subplot(2,1,2)
               plot(record.t_vec,record.dis(:,i),'r','Linewidth',1);
               legendInfo{1} = ('record');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('d(t) [m]');
               xlim([0,40]);
               ylim([-ymax_dis,ymax_dis])
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
                % Comparison plot
               figure(i_fig*i+20+5)
               plot(num_sim.t_vec,num_sim.dis_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               hold on 
               plot(record.t_vec,record.dis(:,i),'r','Linewidth',1);
               legendInfo{2} = ('record');
               grid on
               hold on 
               xlabel('t [s]');
               ylabel('d(t) [m]');
               ylim([-ymax_dis,ymax_dis]);
               xlim([0,40]);
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               legend(legendInfo);
           end
           if flag_RS 
               %% **Comparison: RESPONSE SPECTRA**
               figure(i_fig*i+30)
               plot(record.PSA_T,record.PSA(:,i),'r','Linewidth',1);
               legendInfo{1} = ('record');
               hold on 
               plot(num_sim.PSA_T,num_sim.PSA_orig(:,i),'b','Linewidth',1);
               legendInfo{2} = ('speed original [no filt]');
               grid on
               hold on 
               xlabel('T [s]');
               ylabel('Sa [m/s^2]');
               hold on
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               hold on
               legend(legendInfo);
           end 
       end  
    end
    else
     if flag
       i_fig = 1; 
       for i = 1:3 % e,n,z components 
           %% **Comparison: TIME HISTORIES**
           if flag_TH
               %% ***Comparison: ACCELERATION TH***
               figure(i_fig*i)
               plot(num_sim.t_vec,num_sim.acc_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               xlabel('t [s]');
               ylabel('a(t) [m/s^2]');
               hold on
               title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
               hold on
               legend(legendInfo);
               %% ***Comparison: VELOCITY TH***
               figure(i_fig*i+10)
               plot(num_sim.t_vec,num_sim.vel_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               hold on 
               xlabel('t [s]');
               ylabel('v(t) [m/s]');
               hold on
               title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
               hold on
               legend(legendInfo);
               %% ***Comparison: DISPLACEMNT TH***
               figure(i_fig*i+20)
               plot(num_sim.t_vec,num_sim.dis_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               hold on 
               xlabel('t [s]');
               ylabel('d(t) [m]');
               hold on
               title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
               hold on
               legend(legendInfo);
           end
           if flag_RS 
               %% **Comparison: RESPONSE SPECTRA**
               figure(i_fig*i+30)
               plot(num_sim.PSA_T,num_sim.PSA_orig(:,i),'b','Linewidth',1);
               legendInfo{1} = ('speed original [no filt]');
               hold on 
               xlabel('T [s]');
               ylabel('Sa [m/s^2]');
               hold on
               title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
               hold on
               legend(legendInfo);
           end
  
    end  
     end
end
