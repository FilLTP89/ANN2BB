global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
% _save path_
% sp = '/media/filippo/Data/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images/';
sp = '/media/user/DATI/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images/';
%% *SET UP*
%
% * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
%
pfg.psa = [0 0 16 16];
xlm.psa = {[0,4]};
ylm.psa = {[0,1200]};
xtk.psa = {0:4};
ytk.psa = {0:400:1200};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [g]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
mrk.psa = {'none'};
tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
utd.psa = 100;
%
% * _FOURIER SPECTRA_
%
pfg.fsa = [0 0 16 16];
xlm.fsa = {[1e-2,25]};
ylm.fsa = {[1e-4,1e1]};
xtk.fsa = {[1e-4,1e-3,1e-2,1e-1,1,1e1]};
ytk.fsa = {[1e-2,1e-1,1,10]};
xlb.fsa = {'f [Hz]'};
ylb.fsa = {'FS [m/s]','',''};
scl.fsa = {'log','log','log'};
grd.fsa = {'minor'};
mrk.fsa = {'none'};
tit.fsa = {'FOURIER SPECTRUM'};
utd.fsa = 1;
%
% * _ACCELERATION TIME-HISTORY_
%
pfg.tha.c = [0 0 16 6];
pfg.tha.s = [0 0 16 18];
% xlm.tha = {[0;60];[0;60];[0;60]};
xlm.tha = {[2.5;15];[2.5;15];[2.5;15]};
ylm.tha = {[-1,1];[-1,1];[-1,1]};
% xtk.tha = {0:10:60;0:10:60;0:10:60};
xtk.tha = {[2.5,5:5:15];[2.5,5:5:15];[2.5,5:5:15]};
xlb.tha = {'t [s]'};
ylb.tha = {'a(t) [g]','a(t) [g]','a(t) [g]'};
scl.tha = {'lin','lin','lin'};
grd.tha = {'on'};
mrk.tha = {'none'};
mrk.pga = {'o'};
tit.tha = {'ACCELERATION TIME-HISTORY'};
utd.tha = 1/9.81;
%
% * _VELOCITY TIME-HISTORY_
%
xlm.thv = xlm.tha;
ylm.thv = {[-50,50];[-50,50];[-50,50]};
xtk.thv = xtk.tha;
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
ylm.thd = {[-25,25];[-25,25];[-25,25]};
xtk.thd = xtk.tha;
xlb.thd = xlb.tha;
ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
scl.thd = {'lin','lin','lin'};
grd.thd = {'on'};
mrk.thd = {'none'};
mrk.pgd = {'o'};
tit.thd = {'DISPLACEMENT TIME-HISTORY'};
utd.thd = 100;

cpp = {'e';'n';'z'};

mm_ = 1:4;
syn2ann_plot_res_station;
