%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_drive_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA) to target spectra obtained via
% Artificial Neural Networks.
%% *N.B.*
% Need for:
% _syn2ann_setup.m,syn2ann_rec_drive.m, syn2ann_sim_drive.m,
% syn2ann_emp_sto_drive.m, syn2ann_hybrid_drive.m, syn2ann_ann_drive.m,
% syn2ann_scaling_drive.m, syn2ann_save_res, syn2ann_plot_res_single.m_
%% *REFERENCES*

%% *ANALYSIS SET-UP*
syn2ann_setup;

%% *PARSE RECORDS*
syn2ann_rec_drive;

%% *PARSE NUMERICAL SIMULATIONS*
syn2ann_pbs_drive;

%% *PARSE AND SIMULATE TRAINED ANN*
syn2ann_ann_drive;
fprintf('____________________________________________________________\n');
for NIT=1:MAXIT
    
    fprintf('ITERATION: %u\n',NIT);
    %% *GENERATE EMPIRICAL - PARSE STOCHASTIC*
    syn2ann_emp_sto_drive;
    
    %% *LF-HF CLASSIC HYBRIDIZATION*
    syn2ann_hybrid_drive;
    
    %% *SCORE THE HYBRID TIME HISTORIES*
    if NIT==1
        syn2ann_gof_setup;
    end
    syn2ann_gof_compute;
    
end
fprintf('____________________________________________________________\n');

%% *COMPUTE BEST GOF*
syn2ann_gof_best;

% 
% figure
% plot(pbs.org.mon.vTn,pbs.org.syn{1}.psa.ew,'k');
% hold all;
% for NIT=1:MAXIT
%     plot(pbs.sps{NIT}.mon.vTn(trs.sps.ew.iid),pbs.sps{NIT}.syn{1}.psa.ew(trs.sps.ew.iid),'go--');
%     plot(pbs.sps{NIT}.mon.vTn(trs.sps.ew.tid),pbs.sps{NIT}.syn{1}.psa.ew(trs.sps.ew.tid),'rd--');
% end
% 
% plot(trs.sps.ew.mon.vTn,trs.sps.ew.syn{1}.psa.ew(:),'r');
% 
% figure
% plot(pbs.org.mon.vTn,pbs.org.syn{1}.psa.ew,'k');
% hold all;
% for NIT=1:MAXIT
%     plot(pbs.hyb.sps{NIT}.mon.vTn(trs.sps.ew.iid),pbs.hyb.sps{NIT}.syn{1}.psa.ew(trs.sps.ew.iid),'go--');
%     plot(pbs.hyb.sps{NIT}.mon.vTn(trs.sps.ew.tid),pbs.hyb.sps{NIT}.syn{1}.psa.ew(trs.sps.ew.tid),'rd--');
%     plot(hbs.sps{NIT}.mon.vTn(trs.sps.ew.iid),hbs.sps{NIT}.syn{1}.psa.ew(trs.sps.ew.iid),'cs-');
%     plot(hbs.sps{NIT}.mon.vTn(trs.sps.ew.tid),hbs.sps{NIT}.syn{1}.psa.ew(trs.sps.ew.tid),'bs-');
%     
% end
% 

% col = jet(MAXIT);
% % 
% % 
% % for NIT=1:MAXIT
% %     figure
% %     plot(pbs.org.mon.vTn,pbs.org.syn{1}.psa.ew,'k');
% %     hold all;
% %     plot(trs.sps.ew.mon.vTn,trs.sps.ew.syn{1}.psa.ew(:),'r');
% %     plot(sps.org{NIT}.mon.vTn,sps.org{NIT}.syn{1}.psa.ew,'g');
% %     plot(sps.hyb{NIT}.mon.vTn,sps.hyb{NIT}.syn{1}.psa.ew,'m');
% %     plot(pbs.hyb.sps{NIT}.mon.vTn,pbs.hyb.sps{NIT}.syn{1}.psa.ew,'b');
% %     plot(hbs.sps{NIT}.mon.vTn,hbs.sps{NIT}.syn{1}.psa.ew,'color',col(NIT,:));
% %     keyboard
% % end
% 
% for NIT=1:MAXIT
%     figure
%     loglog(pbs.org.mon.vfr{1},abs(pbs.org.syn{1}.fsa.ew),'k');
%     hold all;
%     loglog(sps.org{NIT}.mon.vfr{1},abs(sps.org{NIT}.syn{1}.fsa.ew),'g');
%     loglog(sps.hyb{NIT}.mon.vfr{1},abs(sps.hyb{NIT}.syn{1}.fsa.ew),'m');
%     keyboard
%     loglog(pbs.hyb.sps{NIT}.mon.vfr{1},abs(pbs.hyb.sps{NIT}.syn{1}.fsa.ew),'b');
%     loglog(hbs.sps{NIT}.mon.vfr{1},abs(hbs.sps{NIT}.syn{1}.fsa.ew),'color',col(NIT,:));
%     
% end
% 
% keyboard



%% *HYB-ANN SPECTRAL MATCHING*
syn2ann_scaling_drive;

%
% %% *HYB-ANN SPECTRAL MATCHING*
% syn2ann_coherency_drive;
%
% %% *SAVE RESULTS*
% syn2ann_save_res;
%
%% *PLOT RESULTS*
syn2ann_plot_res_single;