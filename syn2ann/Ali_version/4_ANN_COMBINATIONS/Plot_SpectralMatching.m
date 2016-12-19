


%% *ORIGINAL HYBRID RESULTS*
t_orig = hybrid.t_vec;
acc_orig = hybrid.acc(:,i);
vel_orig = hybrid.vel(:,i);
dis_orig = hybrid.dis(:,i);
freq_orig = hybrid.freq(:,i);
FAS_orig = hybrid.FAS(:,i);
T_orig = hybrid.PSA_T;
PSA_orig = hybrid.PSA(:,i);
%% *TARGET SPECTRUM*
PSA_tar = PSA_target(:,i);
%% *FINAL BB RESULTS*
t_fin = Out_BB.t_vec;
acc_fin = Out_BB.acc(:,i);
vel_fin = Out_BB.vel(:,i);
dis_fin = Out_BB.dis(:,i);
freq_fin = Out_BB.freq(:,i);
FAS_fin = Out_BB.FAS(:,i);
T_fin = Out_BB.PSA_T;
PSA_fin = Out_BB.PSA(:,i);
%% *COMPARISON PLOT :Time Histories, Fourier Spectra, Response Spectra*
flag = 0;
flag_TH = 1;
flag_FS = 1;
flag_RS = 1;
if flag 
    i_fig = 400;
    %% **Comparison Plot: TIME HISTORIES**
    if flag_TH
        %% ***Comparison Plot: ACCELERATION TH***
        ymax_acc =max(max(acc_orig),max(acc_fin));
        ymin_acc =min(min(acc_orig),min(acc_fin));
        figure(i_fig+i)
        subplot(3,1,1)
%         plot(t_orig,acc_orig,'color',[0.5 0.5 0.5]);
        plot(t_orig,acc_orig,'color','k');
        grid on
        hold on
        xlabel('t [s]'),ylabel('a(t) [m/s^2]');
        xlim([0 45]);
        ylim([ymin_acc,ymax_acc]);
        legend('original hybrid');
        if cfr_record == 1 
           title(strcat(record.station,' - ',record.motion_label(i),' component'));
        else
           title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
        end
        subplot(3,1,2)
        plot(t_fin,acc_fin,'color','blue');
        grid on
        hold on
        xlabel('t [s]');ylabel('a(t) [m/s^2]');
        xlim([0 45]);
        ylim([ymin_acc,ymax_acc]);
        legend('final with ANN');
        subplot(3,1,3)
        plot(t_orig,acc_orig,'color','k');hold on
        plot(t_fin,acc_fin,'color','blue');
        grid on
        xlabel('t [s]');ylabel('a(t) [m/s^2]');
        xlim([0 45]);
        ylim([ymin_acc,ymax_acc]);
        %legend('original hybrid','final with ANN','location','best');
        %% ***Comparison Plot: VELOCITY TH***
        ymax_vel =max(max(vel_orig),max(vel_fin));
        ymin_vel =min(min(vel_orig),min(vel_fin));
        figure(i_fig+10+i)
        subplot(3,1,1)
        plot(t_orig,vel_orig,'color','k');
        grid on
        hold on
        xlabel('t [s]'),ylabel('v(t) [m/s]');
        xlim([0 45]);
        ylim([ymin_vel,ymax_vel]);
        legend('original hybrid');
        if cfr_record == 1 
           title(strcat(record.station,' - ',record.motion_label(i),' component'));
        else
           title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
        end
        subplot(3,1,2)
        plot(t_fin,vel_fin,'color','blue');
        grid on
        hold on
        xlabel('t [s]');ylabel('v(t) [m/s]');
        xlim([0 45]);
        ylim([ymin_vel,ymax_vel]);
        legend('final with ANN');
        subplot(3,1,3)
        plot(t_orig,vel_orig,'color','k');hold on
        plot(t_fin,vel_fin,'color','blue');
        grid on
        xlabel('t [s]');ylabel('v(t) [m/s]');
        xlim([0 45]);
        ylim([ymin_vel,ymax_vel]);
        %legend('original hybrid','final with ANN','location','eastoutside');
        %% ***Comparison Plot: DISPLACEMENT TH***
        figure(i_fig+20+i)
        ymax_dis =max(max(dis_orig),max(dis_fin));
        ymin_dis =min(min(dis_orig),min(dis_fin));
        subplot(3,1,1)
        plot(t_orig,dis_orig,'color','k');
        grid on
        hold on
        xlabel('t [s]'),ylabel('d(t) [m]');
        xlim([0 45]);
        ylim([ymin_dis,ymax_dis]);
        legend('original hybrid');
        if cfr_record == 1 
           title(strcat(record.station,' - ',record.motion_label(i),' component'));
        else
           title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
        end
        subplot(3,1,2)
        plot(t_fin,dis_fin,'color','blue');
        grid on
        hold on
        xlabel('t [s]');ylabel('d(t) [m]');
        xlim([0 45]);
        ylim([ymin_dis,ymax_dis]);
        legend('final with ANN');
        subplot(3,1,3)
        plot(t_orig,dis_orig,'color','k');hold on
        plot(t_fin,dis_fin,'color','blue');
        grid on
        xlabel('t [s]');ylabel('d(t) [m]');
        xlim([0 45]);
        ylim([ymin_dis,ymax_dis]);
        %legend('original hybrid','final with ANN','location','eastoutside');
    end
    %% **Comparison Plot: FOURIER SPECTRA**
    if flag_FS
       figure(i_fig+30+i)
       loglog(freq_orig,FAS_orig,'color','k');hold on
       loglog(freq_fin,FAS_fin,'color','blue');
       grid on
       xlim([0 10]);
       xlabel('f[Hz]');ylabel('FAS [m/s]');
       legend('original','final','location','eastoutside');
       if cfr_record == 1 
          title(strcat(record.station,' - ',record.motion_label(i),' component'));
       else
          title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
       end
    end
    %% **Comparison Plot: RESPONSE SPECTRA**
    if flag_RS
       figure(i_fig+40+i)
       plot(T_orig,PSA_orig,'color','k');hold on
       plot(T_orig,PSA_tar,'color','red'),hold on
       plot(T_orig,(1+tol_upp)*PSA_tar,'--','color','red'),hold on
       plot(T_orig,(1-tol_low)*PSA_tar,'--','color','red'),hold on
       plot(T_fin,PSA_fin,'color','blue');
       legend('original','target_m','target_L_B','target_U_B','final','location','eastoutside');
       xlabel('T [s]'),ylabel('Sa [m/s^2]');
       grid on
       if cfr_record == 1 
          title(strcat(record.station,' - ',record.motion_label(i),' component'));
       else
          title(strcat(num_sim.monID,' - ',num_sim.motion_label(i),' component'));
       end
%     plot(target_Se(:,1),Sp_acc(:),'color',[0.5 0.5 0.5]),hold on
%     plot(target_Se(:,1),target_Se(:,2),'color','red'),hold on
%     plot(target_Se(:,1),(1+tol_upp)*target_Se(:,2),'--','color','red'),hold on
%     plot(target_Se(:,1),(1-tol_low)*target_Se(:,2),'--','color','red'),hold on
%     plot(out_Se(:,1),out_Se(:,2),'color','blue'),hold on
%     legend('original','target_m','target_L_B','target_U_B',...
%         'final','location','eastoutside')
%     xlabel('T(s)'),ylabel('S_e(m/s^2)')
%     
    end
end