%% *ANN COMBINATION*
fprintf('---------------------\n5. FINAL RESULTS\n---------------------\n');

%% *Short description:*
% ...
close all
id_fig = 600;
%% *FINAL COMPARISON PLOT:*
for i = 1:3 % e,n,z components 
    
    T = num_sim.PSA_T;
    % Response Spectrum
    Sa_PBS = num_sim.PSA_orig(:,i);
    Sa_HYB = sel_hyb_PSA(:,i)/100;
    Sa_SPM = Out_BB.PSA(:,i);
    Sa_REC = record.PSA(:,i);
    Sa_targ = PSA_target(:,i);
    % Fourier Amplitude Spectrum
    FAS_PBS = num_sim.FAS_fil(:,i);
    fr_PBS = num_sim.freq_fil(:,i);
    FAS_HYB = hybrid.FAS(:,i);
    fr_HYB = hybrid.freq(:,i);
    FAS_SPM = Out_BB.FAS(:,i);
    fr_SPM = Out_BB.freq(:,i);
    FAS_REC = record.FAS(:,i);
    fr_REC = record.freq(:,i);
    max_fr = [max(fr_PBS);max(fr_REC);max(fr_HYB);max(fr_SPM)];
%     xmax_fr =max(max_fr);
    xmax_fr = 40;
    % Time Histories
    t0=-1;
    t_PBS = num_sim.t_vec+t0;
    acc_PBS = num_sim.acc_orig(:,i);
    vel_PBS = num_sim.vel_orig(:,i);
    dis_PBS = num_sim.dis_orig(:,i);
    t_REC = record.t_vec;
    acc_REC = record.acc(:,i);
    vel_REC = record.vel(:,i);
    dis_REC = record.dis(:,i);
    t_HYB = hybrid.t_vec+t0;
    acc_HYB = sel_hyb_acc(:,i);
    vel_HYB = sel_hyb_vel(:,i);
    dis_HYB = sel_hyb_dis(:,i);
    t_SPM = Out_BB.t_vec+t0;
    acc_SPM = Out_BB.acc(:,i);
    vel_SPM = Out_BB.vel(:,i);
    dis_SPM = Out_BB.dis(:,i);
    
    min_acc = [min(acc_PBS);min(acc_REC);min(acc_HYB);min(acc_SPM)];
    max_acc = [max(acc_PBS);max(acc_REC);max(acc_HYB);max(acc_SPM)];
    min_vel = [min(vel_PBS);min(vel_REC);min(vel_HYB);min(vel_SPM)];
    max_vel = [max(vel_PBS);max(vel_REC);max(vel_HYB);max(vel_SPM)];
    min_dis = [min(dis_PBS);min(dis_REC);min(dis_HYB);min(dis_SPM)];
    max_dis = [max(dis_PBS);max(dis_REC);max(dis_HYB);max(dis_SPM)];
    max_t = [max(t_PBS);max(t_REC);max(t_HYB);max(t_SPM)];
    
    ymin_acc = min(min_acc);
    ymax_acc = max(max_acc);
    ymin_vel = min(min_vel);
    ymax_vel = max(max_vel);
    ymin_dis = min(min_dis);
    ymax_dis = max(max_dis);
    %xmax = max(max_t);
    xmax = 35;
    
    figure(id_fig+i) 
    %% **Final Comparison Plot: RESPONSE SPECTRA**
    subplot(4,4,[4,8]) 
    plot(T,Sa_REC,'r','Linewidth',2);
    hold on 
    plot(T,Sa_PBS,'-.k');
    hold on
    plot(T,Sa_HYB,'--','color',[0.5 0.5 0.5]);
    hold on
    plot(T,Sa_SPM,'k','Linewidth',2);
    hold on
    plot(T,Sa_targ,'g');
    hold on
    legend('REC','PBS','HYB','SPM','Targ');
    xlabel('T [s]'),ylabel('Sa [m/s^2]');
    xlim([0,5]);
    grid on
    if cfr_record == 1 
       title(strcat(record.station,' - ',record.motion_label(i),' component'));
    else
       title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
    end
    %% **Final Comparison Plot: FOURIER SPECTRA**
    subplot(4,4,[12,16])
    loglog(fr_REC,FAS_REC,'r','Linewidth',1);
    hold on 
    loglog(fr_PBS,FAS_PBS,'-.k');
    hold on
    loglog(fr_HYB,FAS_HYB,'color',[0.5 0.5 0.5]);
    hold on
    loglog(fr_SPM,FAS_SPM,'k','Linewidth',1);
    hold on
    xlabel('f [Hz]'),ylabel('FAS [m/s]');
    xlim([0.1,xmax_fr]);
    ylim([10.^(-3),10.^(1)]);
    grid on
    
    %% **Final Comparison Plot: TIME HISTORIES**

    %% ***ACCELERATION TH***
    subplot(4,4,1)
    plot(t_REC,acc_REC,'k');
    hold on
    legend('REC');
    ylim([ymin_acc,ymax_acc]);
    xlim([0,xmax]);
    if cfr_record == 1 
       title(strcat(record.station,' - ',record.motion_label(i),' component'));
    else
       title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
    end
    grid minor
    subplot(4,4,5)
    plot(t_PBS,acc_PBS,'k');
    hold on
    legend('PBS');
    ylim([ymin_acc,ymax_acc]);
    xlim([0,xmax]);
    grid minor
    subplot(4,4,9)
    plot(t_HYB,acc_HYB,'k');
    hold on
    legend('HYB');
    grid minor
    xlim([0,xmax]);
    ylim([ymin_acc,ymax_acc]);
    subplot(4,4,13)
    plot(t_SPM,acc_SPM,'k');
    hold on
    legend('SPM');
    xlim([0,xmax]);
    ylim([ymin_acc,ymax_acc]);
    xlabel('t [s]'),ylabel('a(t) [m/s^2]');
    grid minor
    %% ****VELOCITY TH****
    subplot(4,4,2)
    plot(t_REC,vel_REC,'k');
    hold on
    legend('REC');
    xlim([0,xmax]);
    ylim([ymin_vel,ymax_vel]);
    if cfr_record == 1 
       title(strcat(record.station,' - ',record.motion_label(i),' component'));
    else
       title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
    end
    grid minor
    subplot(4,4,6)
    plot(t_PBS,vel_PBS,'k');
    hold on
    legend('PBS');
    xlim([0,xmax]);
    ylim([ymin_vel,ymax_vel]);
    grid minor
    subplot(4,4,10)
    plot(t_HYB,vel_HYB,'k');
    xlim([0,xmax]);
    ylim([ymin_vel,ymax_vel]);
    hold on
    legend('HYB');
    grid minor
    subplot(4,4,14)
    plot(t_SPM,vel_SPM,'k');
    hold on
    legend('SPM');
    xlabel('t [s]'),ylabel('v(t) [m/s]');
    xlim([0,xmax]);
    ylim([ymin_vel,ymax_vel]);
    grid minor
    %% ****DISPLACEMENT TH****
    subplot(4,4,3)
    plot(t_REC,dis_REC,'k');
    hold on
    legend('REC');
    xlim([0,xmax]);
    ylim([ymin_dis,ymax_dis]);
    if cfr_record == 1 
       title(strcat(record.station,' - ',record.motion_label(i),' component'));
    else
       title(strcat(num_sim.modID,' - ',record.motion_label(i),' component'));
    end
    grid minor
    subplot(4,4,7)
    plot(t_PBS,dis_PBS,'k');
    hold on
    legend('PBS');
    grid minor
    xlim([0,xmax]);
    ylim([ymin_dis,ymax_dis]);
    subplot(4,4,11)
    plot(t_HYB,dis_HYB,'k');
    hold on
    legend('HYB');
    grid minor
    xlim([0,xmax]);
    ylim([ymin_dis,ymax_dis]);
    subplot(4,4,15)
    plot(t_SPM,dis_SPM,'k');
    hold on
    legend('SPM');
    xlabel('t [s]'),ylabel('d(t) [m]');
    xlim([0,xmax]);
    ylim([ymin_dis,ymax_dis]);
    grid minor
    
end 
 