%% *PLOTTING RESULTS*
fprintf('---------------------\n6. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
clr = [rgb('Navy');rgb('IntenseBlue');rgb('Red');rgb('IntenseOrange');rgb('IntenseGreen')];
set(0,'defaultaxescolororder',clr);
cpp = {'e';'n';'z'};
%
% _COMPUTE TIME SHIFT FOR RECORD PLOT_
%
%
% * _FOURIER SPECTRA_
%
pfg.fsa = [0 0 16 16];
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

for mm_ = 1:bhr.ns
    vtm_shift=zeros(1,2);
    switch  bhr.nm{mm_}
        case 'MRN'
            vtm_shift(:) = 1.75*ones(1,2);
            %% *========================= MRN ========================================*
            %
            % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
            %
            pfg.psa = [0 0 16 16];
            xlm.psa = {10.^([log10(0.05),1])};
            ylm.psa = {10.^(1:3:4)};
            xtk.psa = {10.^([log10(0.05),-1:1])};
            ytk.psa = {10.^(1:4)};
            xlb.psa = {'T [s]'};
            ylb.psa = {'PSA [cm/s/s]','',''};
            scl.psa = {'log','log','log'};
            grd.psa = {'minor'};
            mrk.psa = {'none'};
            tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
            utd.psa = 100;
            %
            % * _ACCELERATION TIME-HISTORY_
            %
            pfg.tha.c = [0 0 16 6];
            pfg.tha.s = [0 0 16 18];
            % xlm.tha = {[2;16];[2;16];[2;16]};
            xlm.tha = {[0;40];[0;40];[0;40]};
            ylm.tha = {[-9e2;9e2];[-9e2;9e2];[-9e2;9e2]};
            % xtk.tha = {[2:2:16];[2:2:16];[2:2:16]};
            xtk.tha = {[0:5:40];[0:5:40];[0:5:40]};
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
            ylm.thv = {[-60.0,60.0];[-60.0,60.0];[-60.0,60.0]};
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
            ylm.thd = {[-30.0,30.0];[-30.0,30.0];[-30.0,30.0]};
            xtk.thd = xtk.tha;
            xlb.thd = xlb.tha;
            ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
            scl.thd = {'lin','lin','lin'};
            grd.thd = {'on'};
            mrk.thd = {'none'};
            mrk.pgd = {'o'};
            tit.thd = {'DISPLACEMENT TIME-HISTORY'};
            utd.thd = 100;
        case 'MIR08'
            %% *========================== MIR08 =====================================*
            
            vtm_shift(:) = 0.73*ones(1,2);

            %
            % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
            %
            pfg.psa = [0 0 16 16];
            xlm.psa = {[0,4]};
            ylm.psa = {[0,2000]};
            xtk.psa = {0:4};
            ytk.psa = {0:400:2000};
            xlb.psa = {'T [s]'};
            ylb.psa = {'PSA [cm/s/s]','',''};
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
            % xlm.tha = {[2;16];[2;16];[2;16]};
            xlm.tha = {[0;40];[0;40];[0;40]};
            ylm.tha = {[-5e2,5e2];[-5e2,5e2];[-5e2,5e2]};
            % xtk.tha = {[2:2:16];[2:2:16];[2:2:16]};
            xtk.tha = {[0:5:40];[0:5:40];[0:5:40]};
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
            ylm.thd = {[-18,18];[-18,18];[-18,18]};
            xtk.thd = xtk.tha;
            xlb.thd = xlb.tha;
            ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
            scl.thd = {'lin','lin','lin'};
            grd.thd = {'on'};
            mrk.thd = {'none'};
            mrk.pgd = {'o'};
            tit.thd = {'DISPLACEMENT TIME-HISTORY'};
            utd.thd = 100;
        case 'AQK'
            %% *============================ AQK =====================================*
            
            vtm_shift(:) = 1.25;
            
            %
            % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
            %
            pfg.psa = [0 0 16 16];
            xlm.psa = {[0,4]};
            ylm.psa = {[0,2000]};
            xtk.psa = {0:4};
            ytk.psa = {0:400:2000};
            xlb.psa = {'T [s]'};
            ylb.psa = {'PSA [cm/s/s]','',''};
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
            % xlm.tha = {[0;15];[0;15];[0;15]};
            xlm.tha = {[0;40];[0;40];[0;40]};
            ylm.tha = {[-5e2,5e2];[-5e2,5e2];[-5e2,5e2]};
%             xtk.tha = {[0:3:15];[0:3:15];[0:3:15]};
            xtk.tha = {[0:5:40];[0:5:40];[0:5:40]};
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
            ylm.thd = {[-30,30];[-30,30];[-30,30]};
            xtk.thd = xtk.tha;
            xlb.thd = xlb.tha;
            ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
            scl.thd = {'lin','lin','lin'};
            grd.thd = {'on'};
            mrk.thd = {'none'};
            mrk.pgd = {'o'};
            tit.thd = {'DISPLACEMENT TIME-HISTORY'};
            utd.thd = 100;
        case 'AQU'
            %% *============================ AQU =====================================*
            
            vtm_shift(:) = 1.9*ones(1,2);

            %
            % * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
            %
            pfg.psa = [0 0 16 16];
            xlm.psa = {[0,4]};
            ylm.psa = {[0,2400]};
            xtk.psa = {0:4};
            ytk.psa = {0:400:2400};
            xlb.psa = {'T [s]'};
            ylb.psa = {'PSA [cm/s/s]','',''};
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
            % xlm.tha = {[2;18];[2;18];[2;18]};
            xlm.tha = {[0;40];[0;40];[0;40]};
            ylm.tha = {[-5e2,5e2];[-5e2,5e2];[-5e2,5e2]};
            % xtk.tha = {[2:4:18];[2:4:18];[2:4:18]};
            xtk.tha = {[0:5:40];[0:5:40];[0:5:40]};
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
            ylm.thd = {[-15,15];[-15,15];[-15,15]};
            xtk.thd = xtk.tha;
            xlb.thd = xlb.tha;
            ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
            scl.thd = {'lin','lin','lin'};
            grd.thd = {'on'};
            mrk.thd = {'none'};
            mrk.pgd = {'o'};
            tit.thd = {'DISPLACEMENT TIME-HISTORY'};
            utd.thd = 100;
            %
    end
    syn2ann_plot_res_station;
end