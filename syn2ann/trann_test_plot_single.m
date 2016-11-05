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
xtk.psa = {(0:.5:3)'};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
mrk.psa = {'none'};
tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
utd.psa = 100;

spp = '/home/filippo/Scrivania/ann/';
%spp = fullfile(filesep,'tmp!','gattif','heavy_images_new');

%% *FANCY PLOT*
for mm_ = 1:bhr.ns
    syn2ann_setup_axes;
    trann_fancy_plot;
end
