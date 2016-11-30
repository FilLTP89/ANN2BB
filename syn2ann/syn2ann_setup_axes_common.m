global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

%% *SET-UP*
cpp = bhr.cp;
dlg = cell(3,1);
dlg{1} = 'ew';
dlg{2} = 'ns';
dlg{3} = 'ud';
dlg = upper(intersect(lower(dlg),lower(cpp)));

%% *COLORS*
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

%% *FOURIER SPECTRA*
pfg.fsp = [0 0 28 14];
pfg.fsa = [0 0 10 10];
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

%% * PSEUDO-ACCELERATION RESPONSE SPECTRA
pfg.psa = [0 0 10 10];
xlm.psa = {[0;3]};
xtk.psa = {(0:.5:3)'};
xlb.psa = {'T [s]'};
ylb.psa = {'Sa [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
utd.psa = 100;

%% *TIME-HISTORIES*
pfg.fth = [0 0 28 18.000];
pfg.tha.s = [0 0 14 5.25];
pfg.tha.d = [0 0 14 15.7500];
pfg.tha.t = [0 0 14 47.2500];