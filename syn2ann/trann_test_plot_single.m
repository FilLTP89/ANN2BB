%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_test_plot_single_: function to plot ANN test results at each
% station
%% *N.B.*
% Need for:
% _trann_plot_res_station.m,trann_fancy_plot.m_
fprintf('---------------------\n6. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
% _color sets_
% fancy plot
clr0f = [rgb('Navy');rgb('Red')];
% records-pbs
clr0 = [rgb('Navy');rgb('IntenseBlue')];
% fil-fil-pbs-empirical/stochastic-hybrid
clr10 = [rgb('IntenseBlue');rgb('IntenseGreen');rgb('Orange')];
clr11 = [rgb('DarkGrey');rgb('DarkGrey');rgb('IntenseBlue');rgb('IntenseGreen');rgb('Orange')];
% records-pbs-empirical/stochastic-hybrid
clr12 = [rgb('Navy');rgb('IntenseBlue');rgb('IntenseGreen');rgb('Orange')];
% records-pbs-hybrid-hybrid
clr121 = [rgb('Navy');rgb('Orange');rgb('Magenta')];
% record-hybrid-ann-spectral-matched
clr2 = [rgb('Navy');rgb('Orange');rgb('Red');rgb('Red')];
% records-spectral-matched
clr21 = [rgb('Navy');rgb('Red');rgb('Orange');rgb('IntenseGreen')];
% records-spectral-matched
clr3 = [rgb('Navy');rgb('Red')];
cpp.rec = {'ew';'ns';'ud'};
dlg = cell(3,1);
dlg = {'EW';'NS';'UD'};
pfg.fat = [0 0 30 22];
%
% * _FOURIER SPECTRA_
%
pfg.fsa = [0 0 14 14];
xlm.fsa = {10.^([log10(0.1),log10(40)])};
ylm.fsa = {10.^([-4,1])};
xtk.fsa = {10.^([-1,0,log10(5),1,log10(40)])};
ytk.fsa = {10.^(-4:1)};
xlb.fsa = {'f [Hz]'};
ylb.fsa = {'FS [m/s]','',''};
scl.fsa = {'log','log','log'};
grd.fsa = {'minor'};
mrk.fsa = {'none'};
tit.fsa = {'FOURIER SPECTRUM'};
utd.fsa = 1;
%
% * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
%
pfg.psa = [0 0 14 14];
%     % log-plot
%     xlm.psa = {10.^([log10(0.05);1])};
%     ylm.psa = {10.^(1:3:4)'};
%     xtk.psa = {10.^([log10(0.05);(-1:1)'])};
%     ytk.psa = {10.^(1:4)'};
%     scl.psa = {'log';'log';'log'};
% lin-plot
xlm.psa = {[0;3]};
xtk.psa = {(0:3)'};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
mrk.psa = {'none'};
tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
utd.psa = 100;

spp = '/home/filippo/Scrivania/ann/';

%% *FANCY PLOT*
for mm_ = 1:bhr.ns
    switch  bhr.nm{mm_}
        case 'AMT'
            % *======================== AMT ==============================*
            disp('PLOTTING AMT')
            %             vtm_shift(:) = 0;
            vtm_lim = [5;20];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-6e2;6e2];
            %             thv_lim = [-60;60];
            %             thd_lim = [-30;30];
            ylm.psa = {[0;1200]};
            ytk.psa = {(0:400:1200)'};
        case 'NRC'
            % *======================== NRC ============================*
            disp('PLOTTING NRC')
            %             vtm_shift(:) = 0.73;
            vtm_lim = [5;20];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-4e2;4e2];
            %             thv_lim = [-50;50];
            %             thd_lim = [-20;20];
            ylm.psa = {[0;1800]};
            ytk.psa = {(0:600:1800)'};
        case {'KMM1'}
            % *======================== KMM ==============================*
            disp('PLOTTING KMM')
            %             vtm_shift(:) = 1.25;
            vtm_lim = [15;35];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-4e2;4e2];
            %             thv_lim = [-30;30];
            %             thd_lim = [-30;30];
            ylm.psa = {[0;1200]};
            ytk.psa = {(0:400:1200)'};
        case {'KMM2'}
            % *======================== KMM ==============================*
            disp('PLOTTING KMM')
            %             vtm_shift(:) = 1.25;
            vtm_lim = [12;52];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-10e2;10e2];
            %             thv_lim = [-30;30];
            %             thd_lim = [-30;30];
            ylm.psa = {[0;3000]};
            ytk.psa = {(0:600:3000)'};
    end
    trann_fancy_plot;
end
