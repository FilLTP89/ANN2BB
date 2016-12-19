
%% * COMPARISON: Numerical Simulation vs Synthetic SP96*
flag = 0;
if flag
   i_fig = 100; 
   for i = 1:3 % e,n,z components 
       %% **Comparison: TIME HISTORIES**
       %% ***Comparison: ACCELERATION TH***
       figure(i_fig+i)
       plot(syn_sp96.t_vec,syn_sp96.acc_orig(:,i),'g','Linewidth',1);
       legendInfo{1} = ('sp96 original [no filt]');
       hold on 
       plot(num_sim.t_vec,num_sim.acc_orig(:,i),'b','Linewidth',1);
       legendInfo{2} = ('speed original [no filt]');
       hold on 
       xlabel('t [s]');
       ylabel('a(t) [m/s^2]');
       hold on
       if cfr_record == 1 
          title(strcat(record.station,' - ',record.motion_label(i),' component'));
       else
          title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
       end
       hold on
       legend(legendInfo);
       %% ***Comparison: VELOCITY TH***
       figure(i_fig+i+10)
       plot(syn_sp96.t_vec,syn_sp96.vel_orig(:,i),'g','Linewidth',1);
       legendInfo{1} = ('sp96 original [no filt]');
       hold on 
       plot(num_sim.t_vec,num_sim.vel_orig(:,i),'b','Linewidth',1);
       legendInfo{2} = ('speed original [no filt]');
       hold on 
       xlabel('t [s]');
       ylabel('v(t) [m/s]');
       hold on
       if cfr_record == 1 
          title(strcat(record.station,' - ',record.motion_label(i),' component'));
       else
          title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
       end
       hold on
       legend(legendInfo);
       %% ***Comparison: DISPLACEMNT TH***
       figure(i_fig+i+20)
       plot(syn_sp96.t_vec,syn_sp96.dis_orig(:,i),'g','Linewidth',1);
       legendInfo{1} = ('sp96 original [no filt]');
       hold on 
       plot(num_sim.t_vec,num_sim.dis_orig(:,i),'b','Linewidth',1);
       legendInfo{2} = ('speed original [no filt]');
       hold on 
       xlabel('t [s]');
       ylabel('d(t) [m]');
       hold on
       if cfr_record == 1 
          title(strcat(record.station,' - ',record.motion_label(i),' component'));
       else
          title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
       end
       hold on
       legend(legendInfo);
   end 
end  