global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
%% *SET-UP*

[evt,fgn] = evt2tit(bhr.st{mm_}.ev{1},bhr.st{mm_}.tp{1});
ttt = strcat(bhr.nm{mm_},{' '},evt);
ttt = ttt{1};
evt = evt{1};
st = bhr.nm{mm_};
fprintf('PLOTTING: %s\nDIRECTIONS:\n',ttt);
for i_=1:numel(dlg)
    fprintf('%s\n',dlg{i_});
end
fprintf('\n');
%% *TIME-HISTORIES*
switch  ttt
    case 'MRN 2012-05-29 07:00'
        vtm_shift(:) = 0.5;
        vtm_lim = [0;25];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-9e2;9e2];
        thv_lim = [-60;60];
        thd_lim = [-30;30];
    case 'MIR08'
        vtm_shift(:) = 0.73;
        vtm_lim = [0;25];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-6e2;6e2];
        thv_lim = [-50;50];
        thd_lim = [-20;20];
    case 'AQK-2009-04-06 01:32'
        vtm_shift(:) = 1.35*ones(1,2);
        vtm_lim = [0;25];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        thv_lim = [-30;30];
        thd_lim = [-30;30];
    case 'AQU'
        vtm_shift(:) = 2.15*ones(1,2);
        vtm_lim = [0;25];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        thv_lim = [-30;30];
        thd_lim = [-20;20];
end
%
% * _ACCELERATION TIME-HISTORY_
%

xlm.tha = {vtm_lim;vtm_lim;vtm_lim};
xtk.tha = {vtm_lab;vtm_lab;vtm_lab};
ylm.tha = {tha_lim;tha_lim;tha_lim};
for i_=1:3
    [~,ytk.tha{i_}]=get_axis_tick(xlm.tha{i_},ylm.tha{i_},1,diff(ylm.tha{i_})/4);
end
xlb.tha = {'t [s]'};
ylb.tha = {'a(t) [cm/s/s]','a(t) [cm/s/s]','a(t) [cm/s/s]'};
scl.tha = {'lin','lin','lin'};
grd.tha = {'on'};
mrk.tha = {'none'};
mrk.pga = {'o'};
utd.tha = 100;
%
% * _VELOCITY TIME-HISTORY_
%
xlm.thv = xlm.tha;
ylm.thv = {thv_lim;thv_lim;thv_lim};
xtk.thv = xtk.tha;
for i_=1:3
    [~,ytk.thv{i_}]=get_axis_tick(xlm.thv{i_},ylm.thv{i_},1,diff(ylm.thv{i_})/4);
end
xlb.thv = xlb.tha;
ylb.thv = {'v(t) [cm/s]','v(t) [cm/s]','v(t) [cm/s]'};
scl.thv = {'lin','lin','lin'};
grd.thv = {'on'};
mrk.pgv = {'o'};
utd.thv = 100;
%
% * _DISPLACEMENT TIME-HISTORY_
%
xlm.thd = xlm.tha;
ylm.thd = {thd_lim;thd_lim;thd_lim};
xtk.thd = xtk.tha;
for i_=1:3
    [~,ytk.thd{i_}]=get_axis_tick(xlm.thd{i_},ylm.thd{i_},1,diff(ylm.thd{i_})/4);
end
xlb.thd = xlb.tha;
ylb.thd = {'d(t) [cm]','d(t) [cm]','d(t) [cm]'};
scl.thd = {'lin','lin','lin'};
grd.thd = {'on'};
mrk.thd = {'none'};
mrk.pgd = {'o'};
utd.thd = 100;
%