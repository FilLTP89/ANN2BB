%% *LF-HF HYBRIDIZATION*
if x==1
fprintf('---------------------\n3. LF-HF HYBRIDIZATION\n---------------------\n');
end
%% *Short description:*
% This soubroutine compute the hybrid ground motion using for the LOW
% FREQUENCY (LF) the results of the numerical simulation (num_sim) while
% for the HIGH FREQUENCY (HF) the results of the synthetic obtained with
% Sabetta&Pugliese (SP96). The aim is the computation of a first HYBRID
% version of broadband ground motion (referred as hyb).

%% *RESAMPLING*
%!!!! da valutare se cancellarla forse inutile
dt = num_sim.dt;
for i = 1:3 % e,n,z components 
    [acc_sim_rsmpl(:,i),t_sim_rsmpl,acc_sp96_rsmpl(:,i),t_sp96_rsmpl,dt] = lfhf_rsmpl(num_sim.acc_orig(:,i),num_sim.t_vec,syn_sp96.acc_orig(:,i),syn_sp96.t_vec,dt);
end
%% *PAD LF/HF*
for i = 1:3 % e,n,z components 
    [acc_sim_pad(:,i),acc_sp96_pad(:,i),t_pad] = lfhf_pad(acc_sim_rsmpl(:,i),acc_sp96_rsmpl(:,i),dt);
    vel_sp96_pad = cumtrapz(acc_sp96_pad)*dt;
    dis_sp96_pad = cumtrapz(vel_sp96_pad)*dt;
    vel_sim_pad = cumtrapz(acc_sim_pad)*dt;
    dis_sim_pad = cumtrapz(vel_sim_pad)*dt;
end
%% *ALIGN LF/HF*
for i = 1:3 % e,n,z components 
    t_algn = t_pad;
    [acc_sp96_algn(:,i),vel_sp96_algn(:,i),dis_sp96_algn(:,i)] = lfhf_align(acc_sim_pad(:,i),acc_sp96_pad(:,i),t_pad,dt);
    acc_sim_algn(:,i) = acc_sim_pad(:,i);
    vel_sim_algn(:,i) = vel_sim_pad(:,i);
    dis_sim_algn(:,i) = dis_sim_pad(:,i);
end

check_1 = 1;
if check_1
   for i = 1:3 % e,n,z components 
      % low frequency-arias intensity
      [~,AI5_lf_old(i),~] = arias_intensity(acc_sim_pad(:,i),dt,0.05);
      % high frequency-arias intensity
      [~,AI5_hf_old(i),~] = arias_intensity(acc_sp96_pad(:,i),dt,0.05);
      % low frequency-arias intensity
      [~,AI5_lf_new(i),~] = arias_intensity(acc_sim_algn(:,i),dt,0.05);
      % high frequency-arias intensity
      [~,AI5_hf_new(i),~] = arias_intensity(acc_sp96_algn(:,i),dt,0.05);
      if AI5_lf_new(i)== AI5_hf_new(i)
      else 
         disp ('error');
      end 
  end
end
syn_sp96.acc_orig = zeros(length(acc_sp96_algn(:,1)),3);
syn_sp96.vel_orig = zeros(length(vel_sp96_algn(:,1)),3);
syn_sp96.dis_orig = zeros(length(dis_sp96_algn(:,1)),3);
num_sim.acc_orig = zeros(length(acc_sim_algn(:,1)),3);
num_sim.vel_orig = zeros(length(vel_sim_algn(:,1)),3);
num_sim.dis_orig = zeros(length(dis_sim_algn(:,1)),3);
for i = 1:3
    syn_sp96.acc_orig(:,i) = zeros (size(acc_sp96_algn(:,i)));
    syn_sp96.acc_orig(:,i) = acc_sp96_algn(:,i);
    syn_sp96.vel_orig(:,i) = vel_sp96_algn(:,i);
    syn_sp96.dis_orig(:,i) = dis_sp96_algn(:,i);
    syn_sp96.t_vec = t_algn;
    num_sim.acc_orig(:,i) = acc_sim_algn(:,i);
    num_sim.vel_orig(:,i) = vel_sim_algn(:,i);
    num_sim.dis_orig(:,i) = dis_sim_algn(:,i);
    num_sim.t_vec = t_algn;
end

%% *FILTERING*
%% **Filtering: Define Filtering PARAMETERS*
% Nyquist frequency
fNn = 0.5/dt;
% Cutoff frequency
f_low  = 1.50;
f_high = 1.50;
f_high_up=1.50;
o=4;
%fr=1.5;

%% **Filtering: Filter LOW FREQUENCY records (NUMERICAL SIMULATION SPEED)**
% Butterworth filter's order
bfo = o; 
[bfb,bfa] = butter(bfo,f_low./fNn,'low');
for i = 1:3 % e,n,z components 
    dt = num_sim.dt;
    dis_orig = num_sim.dis_orig(:,i);
    dis_fil = filtfilt(bfb,bfa,dis_orig);
    vel_fil = avd_diff(dt,dis_fil);
    acc_fil = avd_diff(dt,vel_fil);

    num_sim.dis_fil(:,i) = dis_fil;
    num_sim.vel_fil(:,i) = vel_fil;
    num_sim.acc_fil(:,i) = acc_fil;
end
%% **Filtering: Filter HIGH FREQUENCY records (SYNTHETIC SP96)**
% Butterworth filter's order
bfo = o;
[bfb,bfa] = butter(bfo,f_high./fNn,'high');
[bfb_z,bfa_z] = butter(bfo,f_high_up./fNn,'high');
for i = 1:3 % e,n,z components
    dt = syn_sp96.dt;
    dis_orig = syn_sp96.dis_orig(:,i);
    if i==3
    dis_fil = filtfilt(bfb_z,bfa_z,dis_orig);
    else
    dis_fil = filtfilt(bfb,bfa,dis_orig);
    end
    vel_fil = avd_diff(dt,dis_fil);
    acc_fil = avd_diff(dt,vel_fil);
    
    syn_sp96.dis_fil(:,i) = dis_fil;
    syn_sp96.vel_fil(:,i) = vel_fil;
    syn_sp96.acc_fil(:,i) = acc_fil;
end
%% **Filtering: Compute FOURIER TRANSFORM of the Filtered Record LF and HF*
for i = 1:3 % e,n,z components
    %% *** Filtering: FOURIER TRANSFORM of FILTERED NUMERICAL SIMULATION***
    dt = num_sim.dt;
    acc = num_sim.acc_fil(:,i);
    [FT,FAS,freq]= Compute_Fourier(acc,dt);
 
    num_sim.FT_fil(:,i) = FT;
    num_sim.FAS_fil(:,i) = FAS;
    num_sim.freq_fil(:,i) = freq;
    %% *** Filtering: FOURIER TRANSFORM of FILTERED SYNTHETIC SP96***
    dt = syn_sp96.dt;
    acc = syn_sp96.acc_fil(:,i);
    [FT,FAS,freq]= Compute_Fourier(acc,dt);
 
    syn_sp96.FT_fil(:,i) = FT;
    syn_sp96.FAS_fil(:,i) = FAS;
    syn_sp96.freq_fil(:,i) = freq;
end
%% **Filtering: Compute RESPONSE SPECTRA of the Filtered Record LF and HF*
T = num_sim.PSA_T;
syn_sp96.PSA_T = T;
for i = 1:3 % e,n,z components
    for k_=2:length(T)
     %% ***Filtering:RESPONSE SPECTRA of FILTERED NUMERICAL SIMULATION***
        dt = num_sim.dt;
        acc = num_sim.acc_fil(:,i);
        PSA_sim(k_)=newmark_sa(acc,T(k_),0.05,dt);
        PSA_sim(1)=max(abs(acc(:)));
     %% ***Filtering:RESPONSE SPECTRA of FILTERED SYNTHETIC SP96***
        dt = syn_sp96.dt;
        acc = syn_sp96.acc_fil(:,i);
        PSA_syn(k_)=newmark_sa(acc,T(k_),0.05,dt);  
        PSA_syn(1)=max(abs(acc(:)));
    end  
num_sim.PSA_fil(:,i) = PSA_sim;
syn_sp96.PSA_fil(:,i) = PSA_syn;
end

%% *HYBRIDIZATION*
dt = num_sim.dt;
for i = 1:3 % e,n,z components
    %% **Hybridation: Compute HYBRID TIME HISTORIES**
    acc_lf_sim = num_sim.acc_fil(:,i);
    acc_hf_syn = syn_sp96.acc_fil(:,i);
    acc_hb = acc_lf_sim + acc_hf_syn;
    vel_hb = cumtrapz(acc_hb)*dt;
    dis_hb = cumtrapz(vel_hb)*dt;
    %% **Hybridation: Compute HYBRID FOURIER TRANSFORM**
    [FT_hb,FAS_hb,freq_hb]= Compute_Fourier(acc_hb,dt);

    hybrid.dt = num_sim.dt;
    hybrid.t_vec = num_sim.t_vec;
    hybrid.acc(:,i) = acc_hb;
    hybrid.vel(:,i) = vel_hb;
    hybrid.dis(:,i) = dis_hb;
    hybrid.FT(:,i) = FT_hb;
    hybrid.FAS(:,i) = FAS_hb;
    hybrid.freq(:,i) = freq_hb;
end 

%% **Hybridation: Compute RESPONSE SPECTRA**
hybrid.PSA_T = num_sim.PSA_T;
T = hybrid.PSA_T;
dt = hybrid.dt; 

for i = 1:3 % e,n,z components
    acc = hybrid.acc(:,i);
    for k_=2:length(T)
        PSA(k_)=newmark_sa(acc,T(k_),0.05,dt);  
    end  
    PSA(1) = max(abs(acc(:)));
    hybrid.PSA(:,i )= PSA;
end

%% *COMPARISON PLOT: Hybridates, Numerical Simulation, Synthetic sp96, Real Record*
Plot_Hybrid_Comparison;

% %e components
% [acc_sim_e,acc_sp96_e,t_e] = lfhf_pad(acc_sim_e,acc_sp96_e,dt);
% %n components
% [acc_sim_n,acc_sp96_n,t_n] = lfhf_pad(acc_sim_n,acc_sp96_n,dt);
% %z components
% [acc_sim_z,acc_sp96_z,t_z] = lfhf_pad(acc_sim_z,acc_sp96_z,dt);
%!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % low frequency-arias intensity
% [~,AI5_lf_e,~] = arias_intensity(acc_sim_e,dt,0.05);
% % high frequency-arias intensity
% [~,AI5_hf_e,~] = arias_intensity(acc_sp96_e,dt,0.05);
% % low frequency-arias intensity
% [~,AI5_lf_n,~] = arias_intensity(acc_sim_n,dt,0.05);
% % high frequency-arias intensity
% [~,AI5_hf_n,~] = arias_intensity(acc_sp96_n,dt,0.05);
% % low frequency-arias intensity
% [~,AI5_lf_z,~] = arias_intensity(acc_sim_z,dt,0.05);
% % high frequency-arias intensity
% [~,AI5_hf_z,~] = arias_intensity(acc_sp96_z,dt,0.05);
%!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% *ALIGN LF/HF*
% %e components
% [acc_sp96_e,vel_sp96_e,dis_sp96_e] = lfhf_align(acc_sim_e,acc_sp96_e,t_e,dt);
% %n components
% [acc_sp96_n,vel_sp96_n,dis_sp96_n] = lfhf_align(acc_sim_n,acc_sp96_n,t_n,dt);
% %z components
% [acc_sp96_z,vel_sp96_z,dis_sp96_z] = lfhf_align(acc_sim_z,acc_sp96_z,t_z,dt);

%!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % low frequency-arias intensity
% [~,AI5_lf_e,~] = arias_intensity(acc_sim_e,dt,0.05);
% % high frequency-arias intensity
% [~,AI5_hf_e,~] = arias_intensity(acc_sp96_e,dt,0.05);
% % low frequency-arias intensity
% [~,AI5_lf_n,~] = arias_intensity(acc_sim_n,dt,0.05);
% % high frequency-arias intensity
% [~,AI5_hf_n,~] = arias_intensity(acc_sp96_n,dt,0.05);
% % low frequency-arias intensity
% [~,AI5_lf_z,~] = arias_intensity(acc_sim_z,dt,0.05);
% % high frequency-arias intensity
% [~,AI5_hf_z,~] = arias_intensity(acc_sp96_z,dt,0.05);
%!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% *FILTERING*
% % Nnquist frequency
% fNn = 0.5/dt;
% % Cutoff frequency
% fr = 1.5;
% % Butterworth filter's order
% bfo = 2; 
% %% Filter Low Frequencn records (SPEED results)
% [bfb,bfa] = butter(bfo,fr./fNn,'low');
% %% e components [LF records,SPEED results]
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_1 = acc_sim_e;
% % vel_sim_1 = vel_sim_e;
% % dis_sim_1 = dis_sim_e;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dis_sim_e = filtfilt(bfb,bfa,dis_sim_e);
% vel_sim_e = avd_diff(dt,dis_sim_e);
% acc_sim_e = avd_diff(dt,vel_sim_e);
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_2 = acc_sim_e;
% % vel_sim_2 = vel_sim_e;
% % dis_sim_2 = dis_sim_e;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% n components [LF records,SPEED results]
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_1 = acc_sim_n;
% % vel_sim_1 = vel_sim_n;
% % dis_sim_1 = dis_sim_n;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dis_sim_n = filtfilt(bfb,bfa,dis_sim_n);
% vel_sim_n = avd_diff(dt,dis_sim_n);
% acc_sim_n = avd_diff(dt,vel_sim_n);
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_2 = acc_sim_n;
% % vel_sim_2 = vel_sim_n;
% % dis_sim_2 = dis_sim_n;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% z components [LF records,SPEED results]
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_1 = acc_sim_z;
% % vel_sim_1 = vel_sim_z;
% % dis_sim_1 = dis_sim_z;
% %!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dis_sim_z = filtfilt(bfb,bfa,dis_sim_z);
% vel_sim_z = avd_diff(dt,dis_sim_z);
% acc_sim_z = avd_diff(dt,vel_sim_z);
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  acc_sim_2 = acc_sim_z;
% %  vel_sim_2 = vel_sim_z;
% %  dis_sim_2 = dis_sim_z;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Filter High Frequency records (SP96 results)
% [bfb,bfa] = butter(bfo,fr./fNn,'high');
% %% e components [HF records,Sabetta&Pugliese results]
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_1 = acc_sp96_e;
% % vel_sim_1 = vel_sp96_e;
% % dis_sim_1 = dis_sp96_e;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acc_sp96_e = filtfilt(bfb,bfa,acc_sp96_e);
% vel_sp96_e = cumtrapz(acc_sp96_e)*dt;
% dis_sp96_e = cumtrapz(vel_sp96_e)*dt;
% % %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_2 = acc_sp96_e;
% % vel_sim_2 = vel_sp96_e;
% % dis_sim_2 = dis_sp96_e;
% % %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% n components [HF records,Sabetta&Pugliese results]
% %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_1 = acc_sp96_n;
% % vel_sim_1 = vel_sp96_n;
% % dis_sim_1 = dis_sp96_n;
% %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acc_sp96_n = filtfilt(bfb,bfa,acc_sp96_n);
% vel_sp96_n = cumtrapz(acc_sp96_n)*dt;
% dis_sp96_n = cumtrapz(vel_sp96_n)*dt;
% % !!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_2 = acc_sp96_n;
% % vel_sim_2 = vel_sp96_n;
% % dis_sim_2 = dis_sp96_n;
% % !!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% z components [HF records,Sabetta&Pugliese results]
% % %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_1 = acc_sp96_z;
% % vel_sim_1 = vel_sp96_z;
% % dis_sim_1 = dis_sp96_z;
% % %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% acc_sp96_z = filtfilt(bfb,bfa,acc_sp96_z);
% vel_sp96_z = cumtrapz(acc_sp96_z)*dt;
% dis_sp96_z = cumtrapz(vel_sp96_z)*dt;
% % %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % acc_sim_2 = acc_sp96_z;
% % vel_sim_2 = vel_sp96_z;
% % dis_sim_2 = dis_sp96_z;
% % %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% *COMPUTE FOURIER TRANSFORM*
% % Numerical simulation SPEED
% [FT_sim_e,FAS_sim_e,f_sim_e]= Compute_Fourier(acc_sim_e,dt);
% [FT_sim_n,FAS_sim_n,f_sim_n]= Compute_Fourier(acc_sim_n,dt);
% [FT_sim_z,FAS_sim_z,f_sim_z]= Compute_Fourier(acc_sim_z,dt);
% % Synthetic Sabetta&Pugliese '96
% [FT_sp96_e,FAS_sp96_e,f_sp96_e]= Compute_Fourier(acc_sp96_e,dt);
% [FT_sp96_n,FAS_sp96_n,f_sp96_n]= Compute_Fourier(acc_sp96_n,dt);
% [FT_sp96_z,FAS_sp96_z,f_sp96_z]= Compute_Fourier(acc_sp96_z,dt);

% %% *COMPUTE RESPONSE SPECTRA*
% inp_Tn = [0.8:0.1:1.0,1.25:0.25:5.0];
% tar_Tn = [0,0.05,0.1:0.1:0.7,0.75];
% T=[tar_Tn inp_Tn];
% 
% for k_=2:length(T)
%         PSA_sim_e(k_)=newmark_sa(acc_sim_e(:),T(k_),0.05,dt);
%         PSA_sim_n(k_)=newmark_sa(acc_sim_n(:),T(k_),0.05,dt);
%         PSA_sim_z(k_)=newmark_sa(acc_sim_z(:),T(k_),0.05,dt);
%         PSA_sp96_e(k_)=newmark_sa(acc_sp96_e(:),T(k_),0.05,dt);
%         PSA_sp96_n(k_)=newmark_sa(acc_sp96_n(:),T(k_),0.05,dt);
%         PSA_sp96_z(k_)=newmark_sa(acc_sp96_z(:),T(k_),0.05,dt);
% end  
%     PSA_sim_e(1)=max(abs(acc_sim_e(:)));
%     PSA_sim_n(1)=max(abs(acc_sim_n(:)));
%     PSA_sim_z(1)=max(abs(acc_sim_z(:)));
%     PSA_sp96_e(1)=max(abs(acc_sim_e(:)));
%     PSA_sp96_n(1)=max(abs(acc_sim_n(:)));
%     PSA_sp96_z(1)=max(abs(acc_sim_z(:)));


% %% *HBRIDIZATION*
% % e components 
% acc_hb_e = acc_sim_e+acc_sp96_e;
% vel_hb_e = cumtrapz(acc_hb_e)*dt;
% dis_hb_e = cumtrapz(vel_hb_e)*dt;
% [FT_hb_e,FAS_hb_e,f_hb_e,Nfft_hb_e]= Compute_Fourier(acc_hb_e,dt);
% % n components
% acc_hb_n = acc_sim_n+acc_sp96_n;
% vel_hb_n = cumtrapz(acc_hb_n)*dt;
% dis_hb_n = cumtrapz(vel_hb_n)*dt;
% [FT_hb_n,FAS_hb_n,f_hb_n,Nfft_hb_n]= Compute_Fourier(acc_hb_n,dt);
% % z components
% acc_hb_z = acc_sim_z+acc_sp96_z;
% vel_hb_z = cumtrapz(acc_hb_z)*dt;
% dis_hb_z = cumtrapz(vel_hb_z)*dt;
% [FT_hb_z,FAS_hb_z,f_hb_z,Nfft_hb_z]= Compute_Fourier(acc_hb_z,dt);
% % Response spectra e,n,z components
% inp_Tn = [0.8:0.1:1.0,1.25:0.25:5.0];
% tar_Tn = [0,0.05,0.1:0.1:0.7,0.75];
% T=[tar_Tn inp_Tn];
% 
% for k_=2:length(T)
%         PSA_hb_e(k_)=newmark_sa(acc_hb_e(:),T(k_),0.05,dt);
%         PSA_hb_n(k_)=newmark_sa(acc_hb_n(:),T(k_),0.05,dt);
%         PSA_hb_z(k_)=newmark_sa(acc_hb_z(:),T(k_),0.05,dt);
% end  
%     PSA_hb_e(1)=max(abs(acc_hb_e(:)));
%     PSA_hb_n(1)=max(abs(acc_hb_n(:)));
%     PSA_hb_z(1)=max(abs(acc_hb_z(:)));
% t=t_e;

% % %!!!!%%%%%%%%check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag =0;
% if flag
% % e components
% figure(1)
% plot(t_e,acc_hb_e,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('a(m/s^2)');
% hold on 
% plot(t_e,acc_sim_e,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_e,acc_sp96_e,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ew component')
% hold on
% legend(legendInfo);
% 
% figure(11)
% plot(t_e,vel_hb_e,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('v(m/s)');
% hold on 
% plot(t_e,vel_sim_e,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_e,vel_sp96_e,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ew component')
% hold on
% legend(legendInfo);
% 
% figure(12)
% plot(t_e,dis_hb_e,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('d(m)');
% hold on 
% plot(t_e,dis_sim_e,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_e,dis_sp96_e,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ew component')
% hold on
% legend(legendInfo);
% 
% figure(13)
% loglog(f_hb_e,FAS_hb_e,'b');
% legendInfo{1} = ('hybrid');
% xlabel('f(Hz)');
% ylabel('FAS(m/s)');
% hold on 
% plot(f_sim_e,FAS_sim_e,'--r');
% legendInfo{2} = ('speed');
% hold on 
% plot(f_sp96_e,FAS_sp96_e,'--k');
% legendInfo{3} = ['SP96'];
% xlim(10.^([log10(0.05),log10(40)]));
% ylim(10.^([-4,1]));
% hold on 
% % loglog(f_hb_e,[FAS_sim_e+FAS_sp96_e],'g');
% % hold on
% title('MRN ew component')
% hold on
% grid on
% legend(legendInfo);
% 
% figure(14)
% plot(T,PSA_hb_e,'b');
% legendInfo{1} = ('hybrid');
% xlabel('T(s)');
% ylabel('Sa(m/s^2)');
% hold on 
% plot(T,PSA_sim_e,'--r');
% legendInfo{2} = ('speed');
% hold on 
% % xlim(10.^([log10(0.05),log10(40)]));
% % ylim(10.^([-4,1]));
% hold on plot(T,PSA_sp96_e,'--k');
% legendInfo{3} = ('SP96');
% hold on
% title('MRN ew')
% hold on
% grid on
% legend(legendInfo);
% 
% % n components
% figure(2)
% plot(t_n,acc_hb_n,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('a(m/s^2)');
% hold on 
% plot(t_n,acc_sim_n,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_n,acc_sp96_n,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ns component')
% hold on
% legend(legendInfo);
% 
% figure(21)
% plot(t_n,vel_hb_n,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('v(m/s)');
% hold on 
% plot(t_n,vel_sim_n,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_n,vel_sp96_n,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ns component')
% hold on
% legend(legendInfo);
% 
% figure(22)
% plot(t_n,dis_hb_n,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('d(m)');
% hold on 
% plot(t_n,dis_sim_n,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_n,dis_sp96_n,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ns component')
% hold on
% legend(legendInfo);
% 
% figure(23)
% loglog(f_hb_n,FAS_hb_n,'b');
% legendInfo{1} = ('hybrid');
% xlabel('f(Hz)');
% ylabel('FAS(m/s)');
% hold on 
% plot(f_sim_n,FAS_sim_n,'--r');
% legendInfo{2} = ('speed');
% hold on 
% plot(f_sp96_n,FAS_sp96_n,'--k');
% legendInfo{3} = ['SP96'];
% xlim(10.^([log10(0.05),log10(40)]));
% ylim(10.^([-4,1]));
% hold on 
% % loglog(f_hb_n,[FAS_sim_e+FAS_sp96_n],'g');
% % hold on
% title('MRN ns component')
% hold on
% grid on
% legend(legendInfo);
% 
% figure(24)
% plot(T,PSA_hb_n,'b');
% legendInfo{1} = ('hybrid');
% xlabel('T(s)');
% ylabel('Sa(m/s^2)');
% hold on 
% plot(T,PSA_sim_n,'--r');
% legendInfo{2} = ('speed');
% hold on 
% % xlim(10.^([log10(0.05),log10(40)]));
% % ylim(10.^([-4,1]));
% hold on 
% plot(T,PSA_sp96_n,'--k');
% hold on
% legendInfo{3} = ('SP96');
% title('MRN NS')
% hold on
% grid on
% legend(legendInfo);
% 
% % z components
% figure(3)
% plot(t_z,acc_hb_z,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('a(m/s^2)');
% hold on 
% plot(t_z,acc_sim_z,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_z,acc_sp96_z,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ud component')
% hold on
% legend(legendInfo);
% 
% figure(31)
% plot(t_z,vel_hb_z,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('v(m/s)');
% hold on 
% plot(t_z,vel_sim_z,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_z,vel_sp96_z,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ud component')
% hold on
% legend(legendInfo);
% 
% figure(32)
% plot(t_z,dis_hb_z,'b');
% legendInfo{1} = ('hybrid');
% xlabel('t(s)');
% ylabel('d(m)');
% hold on 
% plot(t_z,dis_sim_z,'r');
% legendInfo{2} = ('speed');
% hold on 
% % plot(t_ud,dis_sp96_ud,'--k');
% % legendInfo{3} = ['SP96'];
% % hold on 
% title('MRN ud component')
% hold on
% legend(legendInfo);
% 
% figure(33)
% loglog(f_hb_z,FAS_hb_z,'b');
% legendInfo{1} = ('hybrid');
% xlabel('f(Hz)');
% ylabel('FAS(m/s)');
% hold on 
% plot(f_sim_z,FAS_sim_z,'--r');
% legendInfo{2} = ('speed');
% hold on 
% plot(f_sp96_z,FAS_sp96_z,'--k');
% legendInfo{3} = ['SP96'];
% xlim(10.^([log10(0.05),log10(40)]));
% ylim(10.^([-4,1]));
% hold on 
% % loglog(f_hb_z,[FAS_sim_e+FAS_sp96_z],'g');
% % hold on
% title('MRN ud component')
% hold on
% grid on
% legend(legendInfo);
% 
% figure(34)
% plot(T,PSA_hb_z,'b');
% legendInfo{1} = ('hybrid');
% xlabel('T(s)');
% ylabel('Sa(m/s^2)');
% hold on 
% plot(T,PSA_sim_z,'--r');
% legendInfo{2} = ('speed');
% hold on 
% % xlim(10.^([log10(0.05),log10(40)]));
% % ylim(10.^([-4,1]));
% hold on 
% plot(T,PSA_sp96_z,'--k');
% legendInfo{3} = ('SP96');
% hold on
% title('MRN UD')
% hold on
% grid on
% legend(legendInfo);
% end
% % %!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
