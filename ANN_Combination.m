%% *ANN COMBINATION*
fprintf('---------------------\n4. COMBINATION WITH ANN DATABASE\n---------------------\n');

%% *Short description:*
% This soubroutine computes the broadband ground motion using for the LOW
% FREQUENCY (LF) the results of the numerical simulation (num_sim) while
% for the HIGH FREQUENCY (HF) the results obtained using Artficial Neural
% Network. The aim is the computation of a more accurate version of
% broadband ground motion (correleted to the LF speed results), respect to the first hybrid version.

for i = 1:3 % e,n,z components 
    
    n_=0;

    %% *SET UP PARAMETERS*
    %Tollerance
%     tol_upp = 0.1;
%     tol_low = 0.1;
    tol_upp = 0.5; %Changed by AGO on 15/12
    tol_low = 0.5; %Changed by AGO on 15/12
    %Number of iterations
    niter = 300;   %Changed by AGO on 15/12

    %% *DEFINITION OF THE ORIGINAL HYBRID RESULTS TO BE CORRECTED*
    % Hybrid Acceleration
    hybrid_acc(:,1) = hybrid.t_vec;
    hybrid_acc(:,2) = hybrid.acc(:,i);
    % Hybrid Velocity
    hybrid_vel(:,1) = hybrid.t_vec;
    hybrid_vel(:,2) = hybrid.vel(:,i);
    % Hybrid Displacement
    hybrid_dis(:,1) = hybrid.t_vec;
    hybrid_dis(:,2) = hybrid.dis(:,i);
    % Hybrid Spectra
    hybrid_Se(:,1) = hybrid.PSA_T;
    hybrid_Se(:,2) = hybrid.PSA(:,i);

    %% *DEFINITION OF THE TARGET SPECTRA*
    
    T =[tar_Tn inp_Tn];
    start_inp = find(T==inp_Tn(1));
    end_inp = find(T==inp_Tn(end));
    %PSA(:,i) = hybrid.PSA(:,i);
    PSA(:,i) = num_sim.PSA_orig(:,i);
    PSA_inp(:,i) = PSA(start_inp:end_inp,i);
    PSA_inp(:,i) = PSA_inp(:,i).*100; % cm/s2
    inp(:,i) = log10(PSA_inp(:,i));
    if i == 1||2
         outp(:,i) = sim(net_h,inp(:,i));
         PSA_outp(:,i) = 10.^(outp(:,i));
 %      PSA_outp(:,i) = sim(net_h,PSA_inp(:,i));
    elseif i == 3
        outp(:,i) = sim(net_v,inp(:,i));
        PSA_outp(:,i) = 10.^(outp(:,i));
%       PSA_outp(:,i) = sim(net_v,PSA_inp(:,i));
    end
    PSA_target(:,i) = [PSA_outp(:,i); PSA_inp(:,i)]; % cm/s2
    
    PSA_target(:,i) = PSA_target(:,i);
    target_Se(:,1) = T(:)';
    target_Se(:,2) = PSA_target(:,i);
%     figure(i) 
%     plot(tar_Tn,PSA_outp(:,i),'k');
%     hold on
%     plot([tar_Tn(end),inp_Tn(1)],[PSA_outp(end,i),PSA_inp(1,i)],'--k');
%     hold on
%     grid on
%     plot(inp_Tn,PSA_inp(:,i),'r');
%     ylabel('Sa(cm/s2)');
%     xlabel('T(s)')
%     
 
   %% *SPECTRAL MATCHING*
   % modified by AGO on 15/12
   %[out_acc,out_Se] = SpectralMatching(target_Se,hybrid_acc,hybrid_vel,hybrid_dis,niter,tol_u,tol_l); 
   [out_t,out_acc,out_vel,out_dis,out_T,out_Se,out_freq,out_FAS,n_(i)] = SpectralMatching(target_Se,hybrid_acc,niter,tol_upp,tol_low);
   %[dt,out_acc,out_T,out_Se,out_vel,out_dis] = spectral_matching_Filippo(hybrid.dt,hybrid.acc(:,i),T,PSA_target(:,i));
   %[out_FT,out_FAS,out_freq]= Compute_Fourier(out_acc,dt);
   % end of modification by AGO on 15/12
   %% *CREATE OUTPUT FINAL BB STRUCTURE*
%  Out_BB.t_vec = [0:dt:(size(out_acc,1)-1)*dt]';
   Out_BB.t_vec = out_t;
   Out_BB.acc(:,i) = out_acc./100;
   Out_BB.vel(:,i) = out_vel./100;
   Out_BB.dis(:,i) = out_dis./100;
   Out_BB.PSA_T(:,i) = out_T./100;
   Out_BB.PSA(:,i) = out_Se./100;
   Out_BB.freq(:,i) = out_freq;
   Out_BB.FAS(:,i) = out_FAS./100;
   Out_BB.motion_label(1) = {'e'};
   Out_BB.motion_label(2) = {'n'};
   Out_BB.motion_label(3) = {'z'};
   
   PSA_target(:,i) = PSA_target(:,i)./100;
   
   % modified by AGO on 15/12
   if n_(i)<26
       disp(' ')
       disp('Info: Match is not found to be spectrally satisfactory.')
       disp('      The code will restart, you may want to change the tolerances.') 
       disp('      Press any key to continue...')
       pause
       
       Repeat
   else
       disp('Info: Match is acceptable.')
%        if i==3
%            %% *COMPARISON PLOT: Final Spectral Matched vs Orginal Hybrid*
%            Plot_SpectralMatching;
%    
%            %% *COMPARISON PLOT: Final Spectral Matched vs Real Records*
%            %Plot_FinalSpecMatch_vs_RealRecords
%            clearvars -except wd wd_data wd_results cfr_record record num_sim syn_sp96 hybrid Out_BB PSA_target
%        end
   end

end
   % end of modification by AGO on 15/12    
   


