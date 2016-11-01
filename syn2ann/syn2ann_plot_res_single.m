%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_plot_res_single_: function to plot syn2ann results at each
% station for a single site class
%% *N.B.*
% Need for:
% _get_axis_tick.m,syn2ann_plot_res_station_sp96.m,
% syn2ann_plot_res_station_exsim.m,syn2ann_plot_res_station_both.m,
% syn2ann_fancy_plot.m_
fprintf('---------------------\n7. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
% _color sets_
% fancy plot
clr.rec_pbs_spm = [rgb('Navy');rgb('IntenseBlue');rgb('Red')];
% records-pbs
clr.rec_pbs    = [rgb('Navy');rgb('IntenseBlue')];
% fil-fil-pbs-empirical/stochastic-hybrid
clr.pbs_emp_hyb_ths = [rgb('IntenseBlue');rgb('IntenseGreen');rgb('OrangeRed')];
clr.pbs_emp_hyb_fss = [rgb('DarkGrey');rgb('DarkGrey');rgb('IntenseBlue');rgb('IntenseGreen');rgb('OrangeRed')];
% records-pbs-empirical/stochastic-hybrid
clr.rec_pbs_emp_hyb = [rgb('Navy');rgb('IntenseBlue');rgb('IntenseGreen');rgb('OrangeRed')];
% records-spm-empirical-spm-stochastic-ann
clr.rec_spmemp_spmsto_ann = [rgb('Navy');rgb('IntenseGreen');rgb('Orange');rgb('Red')];
% records-pbs-spm
clr.rec_pbs_spmemp_spmsto = [rgb('Navy');rgb('IntenseBlue');rgb('OrangeRed');rgb('OrangeRed')];
% records-spectral-matched
clr.rec_spm = [rgb('Navy');rgb('Red')];
cpp = {'e';'n';'z'};
pfg.fsp = [0 0 28 14];
pfg.fth = [0 0 28 18.000];
%
% _COMPUTE TIME SHIFT FOR RECORD PLOT_
%
%
% * _FOURIER SPECTRA_
%
pfg.fsa = [0 0 14 14];
xlm.fsa = {10.^([log10(0.1),log10(40)])};
ylm.fsa = {10.^([-3,1])};
xtk.fsa = {10.^([-1,0,log10(5),1,log10(40)])};
ytk.fsa = {10.^(-3:1)};
xlb.fsa = {'f [Hz]'};
ylb.fsa = {'FS [m/s]','',''};
scl.fsa = {'log','log','log'};
grd.fsa = {'minor'};
mrk.fsa = {'none'};
tit.fsa = {'FOURIER SPECTRUM'};
utd.fsa = 1;

for mm_ = 1:bhr.ns
    vtm_shift=zeros(1,2);
    switch  bhr.nm{mm_}
        case 'MRN'
            % *======================== MRN ==============================*
            disp('PLOTTING MRN')
            vtm_shift(:) = 0.5;
            vtm_lim = [0;25];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-9e2;9e2];
            thv_lim = [-60;60];
            thd_lim = [-30;30];
            % lin-plot
            xlm.psa = {[0;3]};
            ylm.psa = {[0;1200]};
            xtk.psa = {(0:3)'};
            ytk.psa = {(0:200:1200)'};
            scl.psa = {'lin','lin','lin'};
            
        case 'MIR08'
            % *======================== MIR08 ============================*
            disp('PLOTTING MIR08')
            vtm_shift(:) = 0.73;
            vtm_lim = [0;25];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-6e2;6e2];
            thv_lim = [-50;50];
            thd_lim = [-20;20];
            % lin-plot
            xlm.psa = {[0;3]};
            ylm.psa = {[0;1200]};
            xtk.psa = {(0:3)'};
            ytk.psa = {(0:200:1200)'};
            scl.psa = {'lin','lin','lin'};
            
        case 'AQK'
            % *======================== AQK ==============================*
            disp('PLOTTING AQK')
            vtm_shift(:) = 1.25;
            vtm_lim = [0;25];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-4e2;4e2];
            thv_lim = [-30;30];
            thd_lim = [-30;30];
            % lin-plot
            xlm.psa = {[0;3]};
            ylm.psa = {[0;1200]};
            xtk.psa = {(0:3)'};
            ytk.psa = {(0:200:1200)'};
            scl.psa = {'lin','lin','lin'};

        case 'AQU'
            % *======================== AQU ==============================*
            disp('PLOTTING AQU')
            vtm_shift(:) = 1.9*ones(1,2);
            vtm_lim = [0;25];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-4e2;4e2];
            thv_lim = [-30;30];
            thd_lim = [-20;20];
    end
    %
    % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
    %
    pfg.psa = [0 0 14 14];
%     %     % log-plot
%     xlm.psa = {10.^([log10(0.05);1])};
%     ylm.psa = {10.^[1,3]};
%     xtk.psa = {10.^[log10(0.05);(-1:1)']};
%     ytk.psa = {10.^(1:3)'};
%     scl.psa = {'log';'log';'log'};
    xlb.psa = {'T [s]'};
    ylb.psa = {'PSA [cm/s/s]','',''};
    grd.psa = {'minor'};
    mrk.psa = {'none'};
    tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
    utd.psa = 100;
    %
    % * _ACCELERATION TIME-HISTORY_
    %
    pfg.tha.c = [0 0 14 5.25];
    pfg.tha.s = [0 0 14 15.7500];
    xlm.tha = {vtm_lim;vtm_lim;vtm_lim};
    xtk.tha = {vtm_lab;vtm_lab;vtm_lab};
    ylm.tha = {tha_lim;tha_lim;tha_lim};
    for iii_=1:3
        [~,ytk.tha{iii_}]=get_axis_tick(xlm.tha{iii_},ylm.tha{iii_},1,diff(ylm.tha{iii_})/4);
    end
    xlb.tha = {'t [s]'};
    ylb.tha = {'a(t) [cm/s/s]','a(t) [cm/s/s]','a(t) [cm/s/s]'};
    scl.tha = {'lin','lin','lin'};
    grd.tha = {'on'};
    mrk.tha = {'none'};
    mrk.pga = {'o'};
    tit.tha = {'ACCELERATION TIME-HISTORY'};
    utd.tha = 100;
    %
    % * _VELOCITY TIME-HISTORY_
    %
    xlm.thv = xlm.tha;
    ylm.thv = {thv_lim;thv_lim;thv_lim};
    xtk.thv = xtk.tha;
    for iii_=1:3
        [~,ytk.thv{iii_}]=get_axis_tick(xlm.thv{iii_},ylm.thv{iii_},1,diff(ylm.thv{iii_})/4);
    end
    xlb.thv = xlb.tha;
    ylb.thv = {'v(t) [cm/s]','v(t) [cm/s]','v(t) [cm/s]'};
    scl.thv = {'lin','lin','lin'};
    grd.thv = {'on'};
    mrk.pgv = {'o'};
    tit.thv = {'VELOCITY TIME-HISTORY'};
    utd.thv = 100;
    %
    % * _DISPLACEMENT TIME-HISTORY_
    %
    xlm.thd = xlm.tha;
    ylm.thd = {thd_lim;thd_lim;thd_lim};
    xtk.thd = xtk.tha;
    for iii_=1:3
        [~,ytk.thd{iii_}]=get_axis_tick(xlm.thd{iii_},ylm.thd{iii_},1,diff(ylm.thd{iii_})/4);
    end
    xlb.thd = xlb.tha;
    ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
    scl.thd = {'lin','lin','lin'};
    grd.thd = {'on'};
    mrk.thd = {'none'};
    mrk.pgd = {'o'};
    tit.thd = {'DISPLACEMENT TIME-HISTORY'};
    utd.thd = 100;
    %
    switch lower(hybrid_type)
        case 'sp96'
            %
            % _SABETTA & PUGLIESE 1996_
            %
            syn2ann_plot_res_station_sp96;
        case 'exsim'
            %
            % _EXSIM_
            %
            syn2ann_plot_res_station_exsim;
        case 'both'
            %
            % _BOTH_
            %
            syn2ann_plot_res_station_both;
    end
%     syn2ann_fancy_plot;
    
end