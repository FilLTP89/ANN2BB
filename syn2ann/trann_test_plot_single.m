%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_plot_res_single_: function to plot syn2ann results at each
% station
%% *N.B.*
% Need for:
% _trann_plot_res_station.m_
fprintf('---------------------\n6. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
% _color sets_
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
dlg{1} = 'EW';
dlg{2} = 'NS';
dlg{3} = 'UD';
%
% _COMPUTE TIME SHIFT FOR RECORD PLOT_
%
%
% * _FOURIER SPECTRA_
%
pfg.fsa = [0 0 14 14];
xlm.fsa = {10.^([log10(0.05),log10(40)])};
ylm.fsa = {10.^([-4,1])};
xtk.fsa = {10.^([log10(0.05),-1,0,log10(5),1,log10(40)])};
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
xlm.psa = {[0.1;5]};
ylm.psa = {[0;2000]};
xtk.psa = {[0.1;(1:5)']};
ytk.psa = {(0:500:2000)'};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
mrk.psa = {'none'};
tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
utd.psa = 100;

vtm_shift=zeros(1,2);
for kk_ = 1:tst.mtd.nr
    cpp.ann = tst.mtd.cpp{kk_};
    flag.ann = seismo_dir_conversion(cpp.ann);
    for mm_ = 1:bhr.ns
        trann_plot_res_station;
    end
end