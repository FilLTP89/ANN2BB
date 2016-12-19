%% *SYNTHETIC SP96*
if x==1
fprintf('---------------------\n2. COMPUTE SYNTHETIC WITH Sabetta&Pugliese (SP96)\n---------------------\n');
end
%% *Short description:*
% This soubroutine compute the synthetic ground motion using the Sabetta&Pugliese 1996, referred as sp96
%% **References:**
% @Article{Sabetta_Pugliese_1996,
%   author =  {Sabetta, F. and Pugliese, A.},
%   title =   {{Estimation of Response Spectra and Simulation of Nonstationary Earthquake Ground Motions}},
%   journal = {Bulletin of the Seismological Society of America},
%   year =    {1996},
%   volume =  {86},
%   number =  {2},
%   pages =   {337--352},

%% *SYNTHETIC RECORD SP96*
% Time step
syn_sp96.dt = num_sim.dt;
%% **Synthetic Record SP96: Compute TIME HISTORIES**
%% ***Synthetic Record SP96: Compute the ACCELERATION TH and time vector***
[t_vec,acc_orig] = sabetta(Mw,R_epi,scc,0,syn_sp96.dt,0.01);
acc_orig = detrend(acc_orig);
syn_sp96.t_vec = t_vec;
for i = 1:3 % e,n,z components
    syn_sp96.acc_orig(:,i) = acc_orig;
end
%% ***Synthetic Record SP96: Compute the VELOCITY and DISPLACEMENT TH***
vel_orig = cumtrapz(acc_orig)*syn_sp96.dt;
dis_orig = cumtrapz(vel_orig)*syn_sp96.dt;
    for i = 1:3 % e,n,z components 
        syn_sp96.vel_orig(:,i) = vel_orig;
        syn_sp96.dis_orig(:,i) = dis_orig;
    end 
syn_sp96.motion_label(1) = {'e'};
syn_sp96.motion_label(2) = {'n'};
syn_sp96.motion_label(3) = {'z'};

%% * COMPARISON: Numerical Simulation vs Synthetic SP96*
Plot_PBS_vs_SP96;

clearvars -except wd wd_data wd_results cfr_record Mw scc R_epi ...
    inp_Tn tar_Tn record num_sim syn_sp96 hybrid net_h net_v ...
    x hyb_EW_acc hyb_NS_acc hyb_Z_acc hyb_EW_PSA hyb_NS_PSA hyb_Z_PSA

% [t_sp96,acc_sp96] = sabetta(Mw,R_epi,scc,0,dt,0.01);
% acc_sp96 = detrend(acc_sp96);
% % velocity
% vel_sp96 = cumtrapz(acc_sp96)*dt;
% % displacement
% dis_sp96 = cumtrapz(vel_sp96)*dt;