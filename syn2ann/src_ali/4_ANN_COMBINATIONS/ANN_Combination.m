%% *ANN COMBINATION*
fprintf('---------------------\n4. COMBINATION WITH ANN DATABASE\n---------------------\n');

%% *Short description:*
% This soubroutine computes the broadband ground motion using for the LOW
% FREQUENCY (LF) the results of the numerical simulation (num_sim) while
% for the HIGH FREQUENCY (HF) the results obtained using Artficial Neural
% Network. The aim is the computation of a more accurate version of
% broadband ground motion (correleted to the LF speed results), 
% respect to the first hybrid version.

% 17/12/16 AGO  : first official release (v2.1)
%                 all unused comments are removed, developers are advised
%                 to refer previous internal alpha versions.

% 18/12/16 AGO  : v3.0 released
%                 code is re-arranged. Different iteration controls are
%                 provided. The automatic loop is removed, since now 500
%                 synthetics are generated.  

% Current version is v. 3.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:3 % e,n,z components 


    %% *SET UP PARAMETERS*
    %Tollerance
    tol_upp_h = 0.15; % Suggested parameter for strict matching
    tol_low_h = 0.15; % Suggested parameter for strict matching
    tol_upp_u = 0.15; % Suggested parameter for strict matching
    tol_low_u = 0.15; % Suggested parameter for strict matching
    tol_upp_pga=0.15; % Suggested parameter for strict matching
    tol_low_pga=0.15; % Suggested parameter for strict matching
    T_sim=0.75;       % Simulation corner spectral period
    
    %Number of iterations
    niter_h = 1000;     % A high number is recommended (Default=300)
    niter_v = 4000;
    
    % loose scheme : 250,  1000
    % medium scheme: 500,  2000
    % stiff scheme : 1000, 4000
    
    %% *DEFINITION OF THE TARGET SPECTRA*
     
    select_the_best;
    
    target_Se(:,1) = T(:)';
    target_Se(:,2) = PSA_target(:,i);
    
   if i<3
       niter=niter_h;
   else
       niter=niter_v;
   end

   %% *DEFINITION OF THE ORIGINAL HYBRID RESULTS TO BE CORRECTED*
    % Hybrid Acceleration
    hybrid_acc(:,1) = hybrid.t_vec;
    hybrid_acc(:,2) = sel_hyb_acc(:,i); %changed here
    % Hybrid Velocity
    hybrid_vel(:,1) = hybrid.t_vec;
    hybrid_vel(:,2) = sel_hyb_vel(:,i); %changed here
    % Hybrid Displacement
    hybrid_dis(:,1) = hybrid.t_vec;
    hybrid_dis(:,2) = sel_hyb_dis(:,i); %changed here
    % Hybrid Spectra
    hybrid_Se(:,1) = hybrid.PSA_T;
    hybrid_Se(:,2) = sel_hyb_PSA(:,i);  %changed here 
    


   %% *SPECTRAL MATCHING*

   if i<2
       tol_upp=tol_upp_h;
       tol_low=tol_low_h;
   else
       tol_upp=tol_upp_u;
       tol_low=tol_low_u;
   end
   
   [out_t,out_acc,out_vel,out_dis,out_T,out_Se,out_freq,out_FAS,n_] = ...
       SpectralMatching(target_Se,hybrid_acc,niter,tol_upp,tol_low,tol_upp_pga,tol_low_pga,i);

   %% *CREATE OUTPUT FINAL BB STRUCTURE*
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
   
   % Checking the convergence criterion
   if n_<27
       disp(' ')
       disp('Info: Match is not found to be spectrally satisfactory.')
       disp('      You may want to observe the agreement, increase the iteration number and re-run.') 
       
       %%% COMMENT IS OPTIONAL, in this version program does not need any
       %%% input from the user.
%        disp('    Press any key to continue...')
%        pause
       %%% COMMENT IS OPTIONAL
       
        %MAIN
   else
       disp('Info: Match is acceptable.')
   end

end
  
   


