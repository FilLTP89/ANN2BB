%% *LOAD INPUT MOTIONS*
fprintf('---------------------\n1. LOAD INPUT MOTIONS\n---------------------\n');

%% *NUMERICAL SIMULATION*

%% **Numerical Simulation: Read INPUT - DISPLACEMENT TH**
% !NON FILTRATI.
% speed_data: time vector; displacement e-w; displacement n-s; displacement ud 
% MARIA ORIGINAL
% if size(num_sim.monID,2) == 5
%    name_file= num_sim.monID;
% elseif size(num_sim.monID,2) == 4
%     name_file= strcat('0',num_sim.monID);
% elseif size(num_sim.monID,2) == 3
%     name_file= strcat('00',num_sim.monID);
% elseif size(num_sim.monID,2) == 2
%     name_file= strcat('000',num_sim.monID);
% elseif size(num_sim.monID,2) == 1
%     name_file= strcat('0000',num_sim.monID);
% end
% MARIA ORIGINAL

%16928-1-FAS-z, 
%18437-2-FAS-ok, 
%12306-FAS-all, 13497, 14120, 15045-done 17869, 18437, 18446
%14120, 18437: matching difficulties, more than one run may needed

speed_data = ...
    importdata('/media/filippo/Data/Filippo/PHD_passing_through_polimi/syn2ann/database/monitors/monitor16928.d');

% time vector
%num_sim.t_vec = speed_data(:,1);
t_vec = speed_data(:,1); %6.12.2016 h. 20.33
% displacement th
for i = 1:3 % e,n,z components 
% num_sim.dis_orig(:,i) = speed_data(:,i+1); %6.12.2016 h. 20.33
dis_orig(:,i) = speed_data(:,i+1);  %6.12.2016 h. 20.33
end
num_sim.motion_label(1) = {'e'};
num_sim.motion_label(2) = {'n'};
num_sim.motion_label(3) = {'z'};
% time step
%num_sim.dt = mean(diff(num_sim.t_vec));%6.12.2016 h. 20.33
num_sim.dt = mean(diff(t_vec));%6.12.2016 h. 20.33
%% **Numerical Simulation: Compute VELOCITY and ACCELERATION TH**
for i = 1:3 % e,n,z components 
    dt = num_sim.dt;
%     dis = num_sim.dis_orig(:,i);
    dis = dis_orig(:,i);
    vel = avd_diff(dt,dis);
    vel = tukeywin(length(vel),0.10).*vel; %6.12.2016 h. 20.33
    vel = [zeros(500,1);vel;zeros(500,1)]; %6.12.2016 h. 20.33
    acc = avd_diff(dt,vel);
    
%     acc = tukeywin(length(acc),0.10).*acc;
%     acc = [zeros(500,1);acc;zeros(500,1)]; %6.12.2016 h. 20.33
%     vel = cumsum(acc).*dt;
%     vel = tukeywin(length(vel),0.10).*vel;
    dis = cumsum(vel).*dt; %6.12.2016 h. 20.33
    t = [0:dt:(length(dis)-1)*dt];%6.12.2016 h. 20.33
    num_sim.t_vec = t;%6.12.2016 h. 20.33
    
    num_sim.dis_orig(:,i) = dis;
    num_sim.vel_orig(:,i) = vel;
    num_sim.acc_orig(:,i) = acc;
end

%% **Numerical Simulation: Compute RESPONSE SPECTRA**
% inp_Tn = [0.8:0.1:1.0,1.25:0.25:5.0]; 
% tar_Tn = [0,0.05,0.1:0.1:0.7,0.75];
T=[tar_Tn inp_Tn];
for i = 1:3 % e,n,z components 
    for k=2:length(T)
        dt = num_sim.dt;
        acc = num_sim.acc_orig(:,i);
        num_sim.PSA_orig(k,i)  = newmark_sa(acc,T(k),0.05,dt);
        num_sim.PSA_T(k)  = T(k);
    end  
    num_sim.PSA_orig(1,i) = max(abs(acc(:)));
    num_sim.PSA_T(1)  = T(1);
end


%% *REAL RECORD*

%% **Real Record: Read INPUT - DISPLACEMENT - VELOCITY - ACCELERATION TH**
if cfr_record == 1
   dis_data = importdata(fullfile(wd_data,'\records\',strcat('record_',record.station,'_dis.d')));
   vel_data = importdata(fullfile(wd_data,'\records\',strcat('record_',record.station,'_vel.d')));
   acc_data = importdata(fullfile(wd_data,'\records\',strcat('record_',record.station,'_acc.d')));
   % time vector
   record.t_vec = dis_data(:,1);
   vtm_shift = 0.5; 
   record.t_vec = record.t_vec + vtm_shift;
   % displacement,velocity,acceleration th
    for i = 1:3 % e,n,z components 
        record.dis(:,i) = dis_data(:,i+1);
        record.vel(:,i) = vel_data(:,i+1);
        record.acc(:,i) = acc_data(:,i+1);
    end
  record.motion_label(1) = {'e'};
  record.motion_label(2) = {'n'};
  record.motion_label(3) = {'z'};
  % time step
  record.dt = mean(diff(record.t_vec));
  %% **Real Record: Compute RESPONSE SPECTRA**
T = num_sim.PSA_T;
for i = 1:3 % e,n,z components 
    for k = 2:length(T)
        dt = record.dt;
        acc = record.acc(:,i);
        record.PSA(k,i)  = newmark_sa(acc,T(k),0.05,dt);
        record.PSA_T(k)  = T(k);
    end  
    record.PSA(1,i) = max(abs(acc(:)));
    record.PSA_T(1)  = T(1);
end
%% **Real Record: Compute FOURIER TRANSFORM**
for i = 1:3 % e,n,z components 
    dt = record.dt;
    acc = record.acc(:,i);
    [FT,FAS,freq]= Compute_Fourier(acc,dt);
    record.FT(:,i) = FT;
    record.FAS(:,i) = FAS;
    record.freq(:,i) = freq;
end
end
%% * COMPARISON: Numerical Simulation vs Real Record*
Plot_PBS_vs_Record;

clearvars -except wd wd_data wd_results cfr_record Mw scc R_epi ...
    inp_Tn tar_Tn record num_sim syn_sp96 hybrid net_h net_v ...
    x hyb_EW_acc hyb_NS_acc hyb_Z_acc hyb_EW_PSA hyb_NS_PSA hyb_Z_PSA
 
% clearvars -except wd wd_data wd_results cfr_record Mw scc R_epi ...
%     inp_Tn tar_Tn record num_sim syn_sp96 hybrid net_h net_v ...
%     x hyb_EW_acc hyb_NS_acc hyb_Z_acc hyb_EW_PSA hyb_NS_PSA hyb_Z_PSA


