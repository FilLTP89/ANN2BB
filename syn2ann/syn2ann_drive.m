%%  *Spectral Matching: Numerical synthetics & ANN*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _syn2ann_drive_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA) to target spectra obtained via
% Artificial Neural Networks.
%% N.B.
% Need for
% _syn2ann_setup.m,syn2ann_rec_drive.m, syn2ann_sim_drive.m,
% syn2ann_sp96_drive.m, syn2ann_exsim_drive.m, syn2ann_hybrid_drive.m,
% syn2ann_ann_drive.m, syn2ann_scaling.m, syn2ann_plot_res.m

%% *SET-UP*
syn2ann_setup;

%% *RECORDS*
syn2ann_rec_drive;
acc = rec.org.syn{1}.tha.e;
dtime=rec.org.mon.dtm(1);
fsa=fft(acc,8192);
close all
figure
loglog(abs(fsa),'b'); hold all;

acc1=detrend(acc);
fsa1=fft(acc1,8192);
loglog(abs(fsa1),'r--'); hold all;

[acc2,~,~]=band_pass_filter(dtime,acc,[],[]);
fsa2=fft(acc2,8192);
loglog(abs(fsa2),'g'); hold all;
keyboard
%% *NUMERICAL SIMULATIONS*
syn2ann_sim_drive;

%% *EMPIRICAL BB SYNTHETICS*
switch lower(hybrid_type)
    case 'sp96'
        % _SABETTA & PUGLIESE 1996_
        syn2ann_sp96_drive;
    case 'exsim'
        % _EXSIM_
        syn2ann_exsim_drive;
end

%% *LF-HF HYBRIDIZATION*
syn2ann_hybrid_drive;

%% *ANN - DATABASE*
syn2ann_ann_drive;

%% *SPECTRAL MATCHING*
syn2ann_scaling_drive;

%% *SAVE RESULTS*
syn2ann_save_res;

%% *PLOT RESULTS*
syn2ann_plot_res_single;