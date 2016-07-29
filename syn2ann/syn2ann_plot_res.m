%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
mm_ = 1;
cpp = {'e';'n';'z'};
%
% _COMPUTE TIME SHIFT FOR RECORD PLOT_
%
syn2ann_alignth;
% vtm_shift(3) = vtm_shift(3)-1.2;
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
%% *========================= MRN ========================================*
mm_ = 1;
%
% * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
%
pfg.psa = [0 0 16 16];
xlm.psa = {[0,4]};
ylm.psa = {[0,2000]};
xtk.psa = {0:4};
ytk.psa = {0:400:2000};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [g]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
mrk.psa = {'none'};
tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
utd.psa = 100;
%
% * _ACCELERATION TIME-HISTORY_
%
pfg.tha.c = [0 0 16 6];
pfg.tha.s = [0 0 16 18];
xlm.tha = {[2;16];[2;16];[2;16]};
ylm.tha = {[-8e2,8e2];[-8e2,8e2];[-8e2,8e2]};
xtk.tha = {[2:2:16];[2:2:16];[2:2:16]};
xlb.tha = {'t [s]'};
ylb.tha = {'a(t) [cm/s^2]','a(t) [cm/s^2]','a(t) [cm/s^2]'};
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
ylm.thv = {[-60,60];[-60,60];[-60,60]};
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

syn2ann_plot_res_station;

% %% *========================== MIR08 =====================================*
% mm_ = 2;
% %
% % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
% %
% pfg.psa = [0 0 16 16];
% xlm.psa = {[0,4]};
% ylm.psa = {[0,1600]};
% xtk.psa = {0:4};
% ytk.psa = {0:400:1600};
% xlb.psa = {'T [s]'};
% ylb.psa = {'PSA [g]','',''};
% scl.psa = {'lin','lin','lin'};
% grd.psa = {'minor'};
% mrk.psa = {'none'};
% tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
% utd.psa = 100;
% %
% % * _ACCELERATION TIME-HISTORY_
% %
% pfg.tha.c = [0 0 16 6];
% pfg.tha.s = [0 0 16 18];
% xlm.tha = {[4;18];[4;18];[4;18]};
% ylm.tha = {[-5e2,5e2];[-5e2,5e2];[-5e2,5e2]};
% xtk.tha = {[4:18];[4:18];[4:18]};
% xlb.tha = {'t [s]'};
% ylb.tha = {'a(t) [cm/s^2]','a(t) [cm/s^2]','a(t) [cm/s^2]'};
% scl.tha = {'lin','lin','lin'};
% grd.tha = {'on'};
% mrk.tha = {'none'};
% mrk.pga = {'o'};
% tit.tha = {'ACCELERATION TIME-HISTORY'};
% utd.tha = 100;
% %
% % * _VELOCITY TIME-HISTORY_
% %
% xlm.thv = xlm.tha;
% ylm.thv = {[-50,50];[-50,50];[-50,50]};
% xtk.thv = xtk.tha;
% xlb.thv = xlb.tha;
% ylb.thv = {'v(t) [cm/s]','v(t) [cm/s]','v(t) [cm/s]'};
% scl.thv = {'lin','lin','lin'};
% grd.thv = {'on'};
% mrk.pgv = {'o'};
% tit.thv = {'VELOCITY TIME-HISTORY'};
% utd.thv = 100;
% %
% % * _DISPLACEMENT TIME-HISTORY_
% %
% xlm.thd = xlm.tha;
% ylm.thd = {[-18,18];[-18,18];[-18,18]};
% xtk.thd = xtk.tha;
% xlb.thd = xlb.tha;
% ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
% scl.thd = {'lin','lin','lin'};
% grd.thd = {'on'};
% mrk.thd = {'none'};
% mrk.pgd = {'o'};
% tit.thd = {'DISPLACEMENT TIME-HISTORY'};
% utd.thd = 100;
% 
% syn2ann_plot_res_station;
% % 
% %% *============================ AQK =====================================*
% mm_ = 3;
% %
% % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
% %
% pfg.psa = [0 0 16 16];
% xlm.psa = {[0,4]};
% ylm.psa = {[0,1600]};
% xtk.psa = {0:4};
% ytk.psa = {0:400:1600};
% xlb.psa = {'T [s]'};
% ylb.psa = {'PSA [g]','',''};
% scl.psa = {'lin','lin','lin'};
% grd.psa = {'minor'};
% mrk.psa = {'none'};
% tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
% utd.psa = 100;
% %
% % * _ACCELERATION TIME-HISTORY_
% %
% pfg.tha.c = [0 0 16 6];
% pfg.tha.s = [0 0 16 18];
% xlm.tha = {[0;15];[0;15];[0;15]};
% ylm.tha = {[-5e2,5e2];[-5e2,5e2];[-5e2,5e2]};
% xtk.tha = {[0:3:15];[0:3:15];[0:3:15]};
% xlb.tha = {'t [s]'};
% ylb.tha = {'a(t) [cm/s^2]','a(t) [cm/s^2]','a(t) [cm/s^2]'};
% scl.tha = {'lin','lin','lin'};
% grd.tha = {'on'};
% mrk.tha = {'none'};
% mrk.pga = {'o'};
% tit.tha = {'ACCELERATION TIME-HISTORY'};
% utd.tha = 100;
% %
% % * _VELOCITY TIME-HISTORY_
% %
% xlm.thv = xlm.tha;
% ylm.thv = {[-60,60];[-60,60];[-60,60]};
% xtk.thv = xtk.tha;
% xlb.thv = xlb.tha;
% ylb.thv = {'v(t) [cm/s]','v(t) [cm/s]','v(t) [cm/s]'};
% scl.thv = {'lin','lin','lin'};
% grd.thv = {'on'};
% mrk.pgv = {'o'};
% tit.thv = {'VELOCITY TIME-HISTORY'};
% utd.thv = 100;
% %
% % * _DISPLACEMENT TIME-HISTORY_
% %
% xlm.thd = xlm.tha;
% ylm.thd = {[-25,25];[-25,25];[-25,25]};
% xtk.thd = xtk.tha;
% xlb.thd = xlb.tha;
% ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
% scl.thd = {'lin','lin','lin'};
% grd.thd = {'on'};
% mrk.thd = {'none'};
% mrk.pgd = {'o'};
% tit.thd = {'DISPLACEMENT TIME-HISTORY'};
% utd.thd = 100;
% 
% syn2ann_plot_res_station;
% % 
% %% *============================ AQU =====================================*
% mm_ = 4;
% %
% % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
% %
% pfg.psa = [0 0 16 16];
% xlm.psa = {[0,4]};
% ylm.psa = {[0,1600]};
% xtk.psa = {0:4};
% ytk.psa = {0:400:1600};
% xlb.psa = {'T [s]'};
% ylb.psa = {'PSA [g]','',''};
% scl.psa = {'lin','lin','lin'};
% grd.psa = {'minor'};
% mrk.psa = {'none'};
% tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
% utd.psa = 100;
% %
% % * _ACCELERATION TIME-HISTORY_
% %
% pfg.tha.c = [0 0 16 6];
% pfg.tha.s = [0 0 16 18];
% xlm.tha = {[2;18];[2;18];[2;18]};
% ylm.tha = {[-5e2,5e2];[-5e2,5e2];[-5e2,5e2]};
% xtk.tha = {[2:4:18];[2:4:18];[2:4:18]};
% xlb.tha = {'t [s]'};
% ylb.tha = {'a(t) [cm/s^2]','a(t) [cm/s^2]','a(t) [cm/s^2]'};
% scl.tha = {'lin','lin','lin'};
% grd.tha = {'on'};
% mrk.tha = {'none'};
% mrk.pga = {'o'};
% tit.tha = {'ACCELERATION TIME-HISTORY'};
% utd.tha = 100;
% %
% % * _VELOCITY TIME-HISTORY_
% %
% xlm.thv = xlm.tha;
% ylm.thv = {[-60,60];[-60,60];[-60,60]};
% xtk.thv = xtk.tha;
% xlb.thv = xlb.tha;
% ylb.thv = {'v(t) [cm/s]','v(t) [cm/s]','v(t) [cm/s]'};
% scl.thv = {'lin','lin','lin'};
% grd.thv = {'on'};
% mrk.pgv = {'o'};
% tit.thv = {'VELOCITY TIME-HISTORY'};
% utd.thv = 100;
% %
% % * _DISPLACEMENT TIME-HISTORY_
% %
% xlm.thd = xlm.tha;
% ylm.thd = {[-15,15];[-15,15];[-15,15]};
% xtk.thd = xtk.tha;
% xlb.thd = xlb.tha;
% ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
% scl.thd = {'lin','lin','lin'};
% grd.thd = {'on'};
% mrk.thd = {'none'};
% mrk.pgd = {'o'};
% tit.thd = {'DISPLACEMENT TIME-HISTORY'};
% utd.thd = 100;
% 
syn2ann_plot_res_station;