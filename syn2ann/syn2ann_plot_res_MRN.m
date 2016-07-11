close all;
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
% _save path_
sp = '/media/filippo/Data/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images/';
% sp = '/media/user/DATI/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images/';
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

nm = 1;

cpp = {'e';'n';'z'};

for i_ = 1:numel(nm)
    st = bhr.nm{nm(i_)};
    lfhf_filter_plot(hbs.mon.dtm(nm(i_)),hbs.mon.vtm{nm(i_)}(end),mon.fa(i_),mon.fb(i_),...
        fullfile(sp,'lfhf_filter'));
    fprintf('%s\n',bhr.nm{nm(i_)});
    for j_ = 1:numel(cpp)
%         %% *POST-PROCESS - ORIGINAL VS FILTERED RECORDS*
%         fn = fullfile(sp,sprintf('%s_rec_org_flt',st));
%         syn2ann_plot_compare(ones(5,1),{rec.org},nm(i_),{'ORIGINAL'},(cpp{j_}),fn,0);
%         
%         %% *POST-PROCESS - NUMERICAL SIMULATIONS VS RECORDS*
%         fn = fullfile(sp,sprintf('%s_sim_rec_org',st));
%         syn2ann_plot_compare(ones(5,1),{rec.org,nss.org},nm(i_),...
%             {'RECORDED';'SIMULATED'},(cpp{j_}),fn,[0,1]);
        
%         %% *POST-PROCESS - BB HYBRIDIZATION*
%         fn = fullfile(sp,sprintf('%s_sim_sp96_hyb',st));
%         syn2ann_plot_compare([0 0 1 1 1],{nss.hyb;sps.hyb;hbs},nm(i_),...
%             {'SIMULATED';'SP96';'HYBRID'},(cpp{j_}),fn,ones(3,1));
%         syn2ann_plot_compare([1 0 0 0 0],{nss.hyb;sps.hyb;hbs},nm(i_),...
%             {'SIMULATED';'SP96';'HYBRID'},(cpp{j_}),fn,ones(3,1));
%         syn2ann_plot_compare([0 1 0 0 0],{nss.org;sps.org;nss.hyb;sps.hyb;hbs},nm(i_),...
%             {'';'';'SIMULATED';'SP96';'HYBRID'},(cpp{j_}),fn,ones(3,1),...
%             {'none','none','none','none','none'},{'--','--','-','-','-'});
        %% *POST-PROCESS - HYBRIDS vs ANN*
        fn = fullfile(sp,sprintf('%s_hyb_ann',st));
        syn2ann_plot_compare([1 0 0 0 0],{hbs;trs.(cpp{j_})},nm(i_),...
            {'HYBRID';'SIMULATED ANN'},(cpp{j_}),fn,[0,0],{'none';'o'},{'-','none'});
        
        %% *POST-PROCESS - SPECTRAL MATCHING vs ANN*
        fn = fullfile(sp,sprintf('%s_spm_ann',st));
        syn2ann_plot_compare([1 0 0 0 0],{spm.(cpp{j_}),trs.(cpp{j_})},nm(i_),...
            {'MATCHED';'SIMULATED ANN'},(cpp{j_}),fn,[0,0],{'none';'o'},{'-','none'});
        
%         %% *POST-PROCESS - SPECTRAL MATCHING vs RECORDS*
%         fn = fullfile(sp,sprintf('%s_spm_rec',st));
%         syn2ann_plot_compare([1 1 1 1 1],{rec.org,spm.(cpp{j_})},nm(i_),...
%             {'RECORDED';'MATCHED'},(cpp{j_}),fn,[0,1]);
        
        %         %% *POST-PROCESS - HYBRIDS VS RECORDS*
        %         fn = fullfile(sp,sprintf('%s_hyb_rec',st));
        %         syn2ann_plot_compare(ones(5,1),{hbs;rec.org},nm(i_),...
        %             {'HYBRID';'RECORDED'},(cpp{j_}),fn,[1,0]);
    end
end