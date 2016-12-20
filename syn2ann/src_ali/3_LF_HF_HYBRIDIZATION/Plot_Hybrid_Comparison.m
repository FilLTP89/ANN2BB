%% *COMPARISON: Hybridates, Numerical Simulation, Synthetic sp96, Real Record*
flag = 0;
flag_TH = 1;
flag_FS = 1;
flag_RS = 1;
   i_fig = 200; 
   for i = 1:3 % e,n,z components
if flag
       if flag_TH
        %% ***Comparison: ACCELERATION TH Hybridates vs Numerical Simulation vs Synthetic SP96**
           ymax_acc = max(hybrid.acc(:,i));
           ymin_acc = min(hybrid.acc(:,i));
           figure(i_fig+i)
           subplot(2,1,1)
           plot(syn_sp96.t_vec,syn_sp96.acc_fil(:,i),'g','Linewidth',1);
           legendInfo{1} = ('sp96 fc 1.5 Hz');
           hold on
           plot(num_sim.t_vec,num_sim.acc_fil(:,i),'m','Linewidth',1);
           legendInfo{2} = ('speed fc 1.5 Hz');
           hold on
           grid on
           legend(legendInfo);
           xlabel('t [s]');
           ylabel('a(t) [m/s^2]');
           xlim([0,40]);
           ylim([ymin_acc,ymax_acc]);
           if cfr_record == 1 
              title(strcat(record.station,' - ',record.motion_label(i),' component'));
           else
              title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
           end
           subplot(2,1,2)
           plot(hybrid.t_vec,hybrid.acc(:,i),'k','Linewidth',1);
           legendInfo{1} = ('hybrid');
           hold on
           legend(legendInfo);
           xlabel('t [s]');
           ylabel('a(t) [m/s^2]');
           xlim([0,40]);
           ylim([ymin_acc,ymax_acc]);
           grid on
           if cfr_record 
            %% ***Comparison: ACCELERATION TH Hybridates vs Real Record**
               ymax_acc = max(max(hybrid.acc(:,i)),max(record.acc(:,i)));
               ymin_acc = min(min(hybrid.acc(:,i)),min(record.acc(:,i)));
               figure(i_fig+i+5)
               subplot(2,1,1)
               plot(hybrid.t_vec,hybrid.acc(:,i),'k','Linewidth',1);
               legendInfo{1} = ('hybrid');
               hold on
               grid on
               xlabel('t [s]');
               ylabel('a(t) [m/s^2]');
               xlim([0,40]);
               ylim([ymin_acc,ymax_acc]);
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               hold on
               legend(legendInfo);
               subplot(2,1,2)
               plot(record.t_vec,record.acc(:,i),'r','Linewidth',1);
               legendInfo{1} = ('record');
               hold on
               grid on
               xlabel('t [s]');
               ylabel('a(t) [m/s^2]');
               xlim([0,40]);
               ylim([ymin_acc,ymax_acc]);
               hold on
               legend(legendInfo);
           end
          %% ***Comparison: VELOCITY TH Hybridates vs Numerical Simulation vs Synthetic SP96**
           ymax_vel = max(hybrid.vel(:,i));
           ymin_vel = min(hybrid.vel(:,i));
           figure(i_fig+i+10)
           subplot(2,1,1)
           plot(syn_sp96.t_vec,syn_sp96.vel_fil(:,i),'g','Linewidth',1);
           legendInfo{1} = ('sp96 fc 1.5 Hz');
           hold on
           plot(num_sim.t_vec,num_sim.vel_fil(:,i),'m','Linewidth',1);
           legendInfo{2} = ('speed fc 1.5 Hz');
           hold on
           grid on
           legend(legendInfo);
           xlabel('t [s]');
           ylabel('v(t) [m/s]');
           xlim([0,40]);
           ylim([ymin_vel,ymax_vel]);
           if cfr_record == 1 
              title(strcat(record.station,' - ',record.motion_label(i),' component'));
           else
              title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
           end
           subplot(2,1,2)
           plot(hybrid.t_vec,hybrid.vel(:,i),'k','Linewidth',1);
           legendInfo{1} = ('hybrid');
           hold on
           legend(legendInfo);
           xlabel('t [s]');
           ylabel('v(t) [m/s]');
           xlim([0,40]);
           ylim([ymin_vel,ymax_vel]);
           grid on
           %% ***Comparison: VELOCITY TH Hybridates vs Real Record**
           ymax_vel = max(max(hybrid.vel(:,i)),max(record.vel(:,i)));
           ymin_vel = min(min(hybrid.vel(:,i)),min(record.vel(:,i)));
           if cfr_record
               figure(i_fig+i+10+5)
               subplot(2,1,1)
               plot(hybrid.t_vec,hybrid.vel(:,i),'k','Linewidth',1);
               legendInfo{1} = ('hybrid');
               hold on
               grid on
               xlabel('t [s]');
               ylabel('v(t) [m/s]');
               xlim([0,40]);
               ylim([ymin_vel,ymax_vel]);
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               hold on
               legend(legendInfo);
               subplot(2,1,2)
               plot(record.t_vec,record.vel(:,i),'r','Linewidth',1);
               legendInfo{1} = ('record');
               hold on
               grid on
               xlabel('t [s]');
               ylabel('v(t) [m/s]');
               xlim([0,40]);
               ylim([ymin_vel,ymax_vel]);
               hold on
               legend(legendInfo);
               hold on
               grid on
               legend(legendInfo);
           end
           %% ***Comparison: DISPLACEMENT TH Hybridates vs Numerical Simulation vs Synthetic SP96**
           ymax_dis = max(hybrid.dis(:,i));
           ymin_dis = min(hybrid.dis(:,i));
           figure(i_fig+i+20)
           subplot(2,1,1)
           plot(num_sim.t_vec,num_sim.dis_fil(:,i),'m','Linewidth',1);
           legendInfo{1} = ('speed fc 1.5 Hz');
           hold on
           plot(syn_sp96.t_vec,syn_sp96.dis_fil(:,i),'g','Linewidth',1);
           legendInfo{2} = ('sp96 fc 1.5 Hz');
           hold on
           grid on
           legend(legendInfo);
           xlabel('t [s]');
           ylabel('d(t) [m]');
           xlim([0,40]);
           ylim([ymin_dis,ymax_dis]);
           if cfr_record == 1 
              title(strcat(record.station,' - ',record.motion_label(i),' component'));
           else
              title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
           end
           subplot(2,1,2)
           plot(hybrid.t_vec,hybrid.dis(:,i),'k','Linewidth',1);
           legendInfo{1} = ('hybrid');
           hold on
           legend(legendInfo);
           xlabel('t [s]');
           ylabel('d(t) [m]');
           xlim([0,40]);
           ylim([ymin_dis,ymax_dis]);
           grid on
           %% ***Comparison: DISPLACEMENT TH Hybridates vs Real Record**
           ymax_dis = max(max(hybrid.dis(:,i)),max(record.dis(:,i)));
           ymin_dis = min(min(hybrid.dis(:,i)),min(record.dis(:,i)));
           if cfr_record
               figure(i_fig+i+20+5)
               plot(hybrid.t_vec,hybrid.dis(:,i),'k','Linewidth',1);
               legendInfo{1} = ('hybrid');
               hold on
               plot(record.t_vec,record.dis(:,i),'r','Linewidth',1);
               legendInfo{2} = ('record');
               grid on
               xlabel('t [s]');
               ylabel('d(t) [m]');
               xlim([0,40]);
               ylim([ymin_dis,ymax_dis]);
               title(strcat(record.station,' - ',record.motion_label(i),' component'));
               hold on
               legend(legendInfo);
            end
       end
       
       if flag_FS
            %% ***Comparison: FOURIER SPECTRA Hybridates vs Numerical Simulation vs Synthetic SP96**
           figure(i_fig+i+30)
           loglog(num_sim.freq_fil,num_sim.FAS_fil(:,i),'m','Linewidth',1);
           legendInfo{1} = ('speed fc 1.5 Hz');
           hold on
           plot(syn_sp96.freq_fil,syn_sp96.FAS_fil(:,i),'g','Linewidth',1);
           legendInfo{2} = ('sp96 fc 1.5 Hz');
           hold on
           plot(hybrid.freq,hybrid.FAS(:,i),'--k','Linewidth',1);
           legendInfo{3} = ('hybrid');
           hold on
           xlabel('f [Hz]');
           ylabel('FAS [m/s]');
           %legend(legendInfo);
           xlim(10.^([log10(0.05),log10(40)]));
           ylim(10.^([-4,1]));
           if cfr_record == 1 
              title(strcat(record.station,' - ',record.motion_label(i),' component'));
           else
              title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
           end
           grid on
           if cfr_record
                %% ***Comparison: FOURIER SPECTRA Hybridates vs Real Record**
               figure(i_fig+i+35)
               loglog(hybrid.freq,hybrid.FAS(:,i),'k','Linewidth',1);
               legendInfo{1} = ('hybrid');
               hold on
               loglog(record.freq,record.FAS(:,i),'r','Linewidth',1);
               legendInfo{2} = ('record');
               hold on
               xlabel('f [Hz]');
               ylabel('FAS [m/s]');
               xlim(10.^([log10(0.05),log10(40)]));
               ylim(10.^([-4,1]));
               if cfr_record == 1 
                  title(strcat(record.station,' - ',record.motion_label(i),' component'));
               else
                  title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
               end
               %legend(legendInfo);
               grid on
           end
           
       end
       %% ***Comparison: RESPONSE SPECTRA Hybridates vs Numerical Simulation vs Synthetic SP96**
       if flag_RS
           figure(i_fig+i+40)
           plot(num_sim.PSA_T,num_sim.PSA_fil(:,i),'m','Linewidth',1);
           legendInfo{1} = ('speed fc 1.5 Hz');
           hold on
           plot(syn_sp96.PSA_T,syn_sp96.PSA_fil(:,i),'g','Linewidth',1);
           legendInfo{2} = ('sp96 fc 1.5 Hz');
           hold on
           plot(hybrid.PSA_T,hybrid.PSA(:,i),'k','Linewidth',1);
           legendInfo{3} = ('hybrid');
           hold on
           xlabel('T [s]');
           ylabel('Sa [m/s^2]');
           if cfr_record == 1 
              title(strcat(record.station,' - ',record.motion_label(i),' component'));
           else
              title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
           end
           hold on
           legend(legendInfo);
           grid on
        %% ***Comparison: RESPONSE SPECTRA Hybridates vs Real Record**
        if cfr_record
           figure(i_fig+i+45)
           plot(hybrid.PSA_T,hybrid.PSA(:,i),'k','Linewidth',1);
           legendInfo{1} = ('hybrid');
           hold on 
           plot(record.PSA_T,record.PSA(:,i),'r','Linewidth',1);
           legendInfo{2} = ('real record');
           hold on
           xlabel('T [s]');
           ylabel('Sa [m/s^2]');
           if cfr_record == 1 
              title(strcat(record.station,' - ',record.motion_label(i),' component'));
           else
              title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
           end
           hold on
           legend(legendInfo);
           grid on
        end
  end
end

end
