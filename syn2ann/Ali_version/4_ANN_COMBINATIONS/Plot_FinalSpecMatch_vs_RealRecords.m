
if cfr_record == 1
    %% *FINAL BB RESULTS*
    t_fin = Out_BB.t_vec;
    acc_fin = Out_BB.acc(:,i);
    vel_fin = Out_BB.vel(:,i);
    dis_fin = Out_BB.dis(:,i);
    freq_fin = Out_BB.freq(:,i);
    FAS_fin = Out_BB.FAS(:,i);
    T_fin = Out_BB.PSA_T;
    PSA_fin = Out_BB.PSA(:,i);

    %% * REAL RECORDS*
    t_rec = record.t_vec;
    acc_rec = record.acc(:,i);
    vel_rec = record.vel(:,i);
    dis_rec = record.dis(:,i);
    freq_rec = record.freq(:,i);
    FAS_rec = record.FAS(:,i);
    T_rec = record.PSA_T;
    PSA_rec = record.PSA(:,i);
    
    flag = 0;
    flag_TH = 1;
    flag_FS = 1;
    flag_RS = 1;
    
if flag
   i_fig = 500;
   %% **Comparison Plot: TIME HISTORIES**
   if flag_TH
      %% ***Comparison Plot: ACCELERATION TH***
      ymax_acc=max(max(acc_fin),max(acc_rec));
      ymin_acc=min(min(acc_fin),min(acc_rec));
      figure(i_fig+i)
      subplot(2,1,1)
      plot(t_fin,acc_fin,'color','blue');
      grid on
      hold on
      xlabel('t [s]'),ylabel('a(t) [m/s^2]');
      xlim([0 45]);
      ylim([ymin_acc,ymax_acc]);
      legend('final simulated BB');
      if cfr_record == 1 
         title(strcat(record.station,' - ',record.motion_label(i),' component'));
      else
         title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
      end
      subplot(2,1,2)
      plot(t_rec,acc_rec,'color','red');
      xlim([0 45]);
      grid on
      hold on
      xlabel('t [s]');ylabel('a(t) [m/s^2]');
      xlim([0 45]);
      ylim([ymin_acc,ymax_acc]);
      legend('real record');
      %% ***Comparison Plot: VELOCITY TH***
      ymax_vel=max(max(vel_fin),max(vel_rec));
      ymin_vel=min(min(vel_fin),min(vel_rec));
      figure(i_fig+10+i)
      subplot(2,1,1)
      plot(t_fin,vel_fin,'color','blue');
      grid on
      hold on
      xlabel('t [s]'),ylabel('v(t) [m/s]');
      xlim([0 45]);
      ylim([ymin_vel,ymax_vel]);
      legend('final simulated BB');
      if cfr_record == 1 
         title(strcat(record.station,' - ',record.motion_label(i),' component'));
      else
         title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
      end
      subplot(2,1,2)
      plot(t_rec,vel_rec,'color','red');
      grid on
      hold on
      xlabel('t [s]');ylabel('v(t) [m/s]');
      xlim([0 45]);
      ylim([ymin_vel,ymax_vel]);
      legend('real record');
      %% ***Comparison Plot: DISPLACEMENT TH***
      ymax_dis=max(max(dis_fin),max(dis_rec));
      ymin_dis=min(min(dis_fin),min(dis_rec));
      figure(i_fig+20+i)
      plot(t_fin,dis_fin,'color','blue');
      hold on 
      plot(t_rec,dis_rec,'color','red');
      grid on
      hold on
      xlabel('t [s]'),ylabel('d(t) [m]');
      xlim([0 45]);
      ylim([ymin_dis,ymax_dis]);
      legend('final simulated BB','real record');
      if cfr_record == 1 
         title(strcat(record.station,' - ',record.motion_label(i),' component'));
      else
         title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
      end
   end
  if flag_FS
     figure(i_fig+30+i)
     loglog(freq_rec,FAS_rec,'r');hold on
     loglog(freq_fin,FAS_fin,'color','blue');
     grid on
     xlim([0 10]);
     xlabel('f[Hz]');ylabel('FAS [m/s]');
     legend('real record','final simulated BB','location','eastoutside');
     if cfr_record == 1 
        title(strcat(record.station,' - ',record.motion_label(i),' component'));
     else
        title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
     end
    end
    if flag_RS
       figure(i_fig+40+i)
       plot(T_rec,PSA_rec,'r');hold on
       plot(T_fin,PSA_fin,'color','blue');
       legend('real record','final simulated BB','location','eastoutside');
       xlabel('T [s]'),ylabel('Sa [m/s^2]');
       grid on
       if cfr_record == 1 
          title(strcat(record.station,' - ',record.motion_label(i),' component'));
       else
          title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
       end
    end
    end
    
end