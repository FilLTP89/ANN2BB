%% *PLOTTING RESULTS*
fprintf('---------------------\n6. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
clr = [rgb('Navy');rgb('Red');rgb('IntenseGreen')];
clr1 = [rgb('Navy');rgb('Red');rgb('Grey');rgb('Grey');rgb('IntenseGreen')];
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
            disp('PLOTTING MRN')
            
            %% *========================= MRN ========================================*
            vtm_shift(:) = 1.75*ones(1,2);
            
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
            disp('PLOTTING MIR08')
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
            disp('PLOTTING AQK')
            %% *============================ AQK =====================================*
            
            vtm_shift(:) = 1.25;
            vtm_lim = [0;25];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-4e2;4e2];
            thv_lim = [-30;30];
            thd_lim = [-30;30];
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
            xlm.tha = {vtm_lim;vtm_lim;vtm_lim};
            ylm.tha = {tha_lim;tha_lim;tha_lim};
            xtk.tha = {vtm_lab;vtm_lab;vtm_lab};
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
            xlb.thd = xlb.tha;
            ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
            scl.thd = {'lin','lin','lin'};
            grd.thd = {'on'};
            mrk.thd = {'none'};
            mrk.pgd = {'o'};
            tit.thd = {'DISPLACEMENT TIME-HISTORY'};
            utd.thd = 100;
            
        case 'AQU'
            disp('PLOTTING AQU')
            %% *============================ AQU =====================================*
            
            vtm_shift(:) = 1.9*ones(1,2);
            vtm_lim = [0;25];
            vtm_lab = (vtm_lim(1):5:vtm_lim(end));
            tha_lim = [-4e2;4e2];
            thv_lim = [-30;30];
            thd_lim = [-20;20];
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
            xlm.tha = {vtm_lim;vtm_lim;vtm_lim};
            ylm.tha = {tha_lim;tha_lim;tha_lim};
            xtk.tha = {vtm_lab;vtm_lab;vtm_lab};
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
            xlb.thd = xlb.tha;
            ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
            scl.thd = {'lin','lin','lin'};
            grd.thd = {'on'};
            mrk.thd = {'none'};
            mrk.pgd = {'o'};
            tit.thd = {'DISPLACEMENT TIME-HISTORY'};
            utd.thd = 100;
            keyboard
            %
    end
    
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
    end
    
end