close all;
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit
% _save path_
% sp = '/media/filippo/Data/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images/';
sp = '/media/user/DATI/Filippo/PHD_passing_through/WRITEUP/ALL_PICTURES/images/';
%% *SET UP*
%
% * _PSEUDO-ACCELERATION RESPONSE SPECTRA_
%
pfg.psa = [0 0 16 16];
xlm.psa = {[0,5]};
ylm.psa = {[0,25]};
xtk.psa = {[0:5]};
ytk.psa = {[0:5:25]};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [m/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
mrk.psa = {'none'};
tit.psa = {'PSEUDO-ACCELERATION SPECTRUM (5%)'};
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
%
% * _ACCELERATION TIME-HISTORY_
%
pfg.tha.c = [0 0 16 6];
pfg.tha.s = [0 0 16 18];
xlm.tha = {[0;30];[0;30];[0;30]};
xtk.tha = {[0:5:30];[0:5:30];[0:5:30]};
xlb.tha = {'','','t [s]'};
ylb.tha = {'a(t) [m/s/s]','a(t) [m/s/s]','a(t) [m/s/s]'};
scl.tha = {'lin','lin','lin'};
grd.tha = {'on'};
mrk.tha = {'none'};
mrk.pga = {'o'};
tit.tha = {'ACCELERATION TIME-HISTORY'};
%
% * _VELOCITY TIME-HISTORY_
%
xlm.thv = {[0;20];[0;20];[0;20]};
xtk.thv = {[0:5:20];[0:5:20];[0:5:20]};
xlb.thv = {'','','t [s]'};
ylb.thv = {'v(t) [m/s]','v(t) [m/s]','v(t) [m/s]'};
scl.thv = {'lin','lin','lin'};
grd.thv = {'on'};
mrk.pgv = {'o'};
tit.thv = {'VELOCITY TIME-HISTORY'};
%
% * _DISPLACEMENT TIME-HISTORY_
%
xlm.thd = {[0;20];[0;20];[0;20]};
xtk.thd = {[0:5:20];[0:5:20];[0:5:20]};
xlb.thd = {'','','t [s]'};
ylb.thd = {'d(t) [m]','d(t) [m]','d(t) [m]'};
scl.thd = {'lin','lin','lin'};
grd.thd = {'on'};
mrk.thd = {'none'};
mrk.pgd = {'o'};
tit.thd = {'DISPLACEMENT TIME-HISTORY'};

nm = 1;
st  = strcat(bhr.st{1}.id,bhr.st{1}.dv{1}); st = st{1};
cpp = {'n','z'};

for i_ = 1:numel(nm)
    for j_ = 1:numel(cpp)
        %% *POST-PROCESS - ORIGINAL VS FILTERED RECORDS*
        fn = fullfile(sp,sprintf('%s_rec_org_flt',st));
        syn2ann_plot_compare(ones(5,1),{rec.org},nm(i_),{'ORIGINAL';'FILTERED'},(cpp{j_}),fn,0);
        %% *POST-PROCESS - NUMERICAL SIMULATIONS VS RECORDS*
        fn = fullfile(sp,sprintf('%s_sim_rec_org',st));
        syn2ann_plot_compare(ones(5,1),{nss.org;rec.org},nm(i_),{'SIMULATED';'RECORDED'},(cpp{j_}),fn,[1,0]);
        %% *POST-PROCESS - BB HYBRIDIZATION*
        fn = fullfile(sp,sprintf('%s_sim_sp96_hyb',st));
        syn2ann_plot_compare(ones(5,1),{nss.hyb;sps.hyb;hbs},nm(i_),{'SIMULATED';'SP96';'HYBRID'},(cpp{j_}),fn,[0,0,0]);
        %% *POST-PROCESS - HYBRIDS VS RECORDS*
        fn = fullfile(sp,sprintf('%s_hyb_rec',st));
        syn2ann_plot_compare(ones(5,1),{hbs;rec.org},nm(i_),{'HYBRID';'RECORDED'},(cpp{j_}),fn,[1,0]);
        %% *POST-PROCESS - HYBRIDS vs ANN*
        fn = fullfile(sp,sprintf('%s_hyb_ann',st));
        syn2ann_plot_compare([1 0 0 0 0],{hbs;trs.(cpp{j_})},nm(i_),{'HYBRID';'SIMULATED ANN'},(cpp{j_}),fn,[0,0]);
        %% *POST-PROCESS - SPECTRAL MATCHING vs RECORDS*
        fn = fullfile(sp,sprintf('%s_spm_rec',st));
        syn2ann_plot_compare([1 0 1 1 1],{spm.(cpp{j_}),rec.org},nm(i_),{'SPECTRAL MATCHED';'RECORDED'},(cpp{j_}),fn,[1,0]);
        %% *POST-PROCESS - SPECTRAL MATCHING vs ANN*
        fn = fullfile(sp,sprintf('%s_spm_ann',st));
        syn2ann_plot_compare([1 0 0 0 0],{spm.(cpp{j_}),trs.(cpp{j_})},nm(i_),...
            {'SPECTRAL MATCHED';'SIMULATED ANN'},(cpp{j_}),fn,[0,0]);
    end
end