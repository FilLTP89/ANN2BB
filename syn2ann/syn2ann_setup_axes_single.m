global xlm xtk ylm ytk utd
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
        xlm.tha = [0;25];
        xtk.tha = (xlm.tha(1):5:xlm.tha(end))';
    case 'MIR08'
        vtm_shift(:) = 0.73;
        xlm.tha = [0;25];
        xtk.tha = (xlm.tha(1):5:xlm.tha(end))';
    case 'AQK 2009-04-06 01:32'
        vtm_shift(:) = 1.35;
        xlm.tha = [0;25];
        xtk.tha = (xlm.tha(1):5:xlm.tha(end))';
    case 'AQU'
        vtm_shift(:) = 2.15;
        xlm.tha = [0;25];
        xtk.tha = (xlm.tha(1):5:xlm.tha(end))';
end
xlm.thv = xlm.tha;
xtk.thv = xtk.tha;
xlm.thd = xlm.tha;
xtk.thd = xtk.tha;

%% *DEFINE TIME-HISTORY YLIM*
tha_lim = 0;
thv_lim = 0;
thd_lim = 0;
psa_lim = 0;
for j_ = 1:numel(cpp)
    tha_lim = max([tha_lim;abs(rec.org.syn{mm_}.pga.(cpp{j_})(2))]);
    thv_lim = max([thv_lim;abs(rec.org.syn{mm_}.pgv.(cpp{j_})(2))]);
    thd_lim = max([thd_lim;abs(rec.org.syn{mm_}.pgd.(cpp{j_})(2))]);
    tha_lim = max([tha_lim;abs(nss.org.syn{mm_}.pga.(cpp{j_})(2))]);
    thv_lim = max([thv_lim;abs(nss.org.syn{mm_}.pgv.(cpp{j_})(2))]);
    thd_lim = max([thd_lim;abs(nss.org.syn{mm_}.pgd.(cpp{j_})(2))]);
    psa_lim = max([psa_lim;max(abs(rec.org.syn{mm_}.psa.(cpp{j_})))]);
end
switch lower(hybrid_type)
    case 'sp96'
        for j_ = 1:numel(cpp)
            tha_lim = max([tha_lim;abs(hbs.sps.syn{mm_}.pga.(cpp{j_})(2))]);
            thv_lim = max([thv_lim;abs(hbs.sps.syn{mm_}.pgv.(cpp{j_})(2))]);
            thd_lim = max([thd_lim;abs(hbs.sps.syn{mm_}.pgd.(cpp{j_})(2))]);
            psa_lim = max([psa_lim;max(abs(hbs.sps.syn{mm_}.psa.(cpp{j_})))]);
            tha_lim = max([tha_lim;abs(spm.sps.(cpp{j_}).syn{mm_}.pga.(cpp{j_})(2))]);
            thv_lim = max([thv_lim;abs(spm.sps.(cpp{j_}).syn{mm_}.pgv.(cpp{j_})(2))]);
            thd_lim = max([thd_lim;abs(spm.sps.(cpp{j_}).syn{mm_}.pgd.(cpp{j_})(2))]);
            psa_lim = max([psa_lim;max(abs(spm.sps.(cpp{j_}).syn{mm_}.psa.(cpp{j_})))]);
        end
    case 'exsim'
        for j_ = 1:numel(cpp)
            tha_lim = max([tha_lim;abs(hbs.exs.syn{mm_}.pga.(cpp{j_})(2))]);
            thv_lim = max([thv_lim;abs(hbs.exs.syn{mm_}.pgv.(cpp{j_})(2))]);
            thd_lim = max([thd_lim;abs(hbs.exs.syn{mm_}.pgd.(cpp{j_})(2))]);
            psa_lim = max([psa_lim;max(abs(hbs.exs.syn{mm_}.psa.(cpp{j_})))]);
            tha_lim = max([tha_lim;abs(spm.exs.(cpp{j_}).syn{mm_}.pga.(cpp{j_})(2))]);
            thv_lim = max([thv_lim;abs(spm.exs.(cpp{j_}).syn{mm_}.pgv.(cpp{j_})(2))]);
            thd_lim = max([thd_lim;abs(spm.exs.(cpp{j_}).syn{mm_}.pgd.(cpp{j_})(2))]);
            psa_lim = max([psa_lim;max(abs(spm.exs.(cpp{j_}).syn{mm_}.psa.(cpp{j_})))]);
        end
end
tha_lim = tha_lim*utd.tha;
thv_lim = thv_lim*utd.thv;
thd_lim = thd_lim*utd.thd;
psa_lim = psa_lim*utd.psa;

ylm.tha = [-ceil(tha_lim/4)*4,ceil(tha_lim/4)*4];
[~,ytk.tha] = get_axis_tick(ylm.tha,ylm.tha,ceil(ylm.tha(2)/2),ceil(ylm.tha(2)/2));
ylm.tha = ylm.tha;
ytk.tha = ytk.tha;

ylm.thv = [-ceil(thv_lim/4)*4,ceil(thv_lim/4)*4];
[~,ytk.thv] = get_axis_tick(ylm.thv,ylm.thv,ceil(ylm.thv(2)/2),ceil(ylm.thv(2)/2));
ylm.thv = ylm.thv;
ytk.thv = ytk.thv;

ylm.thd = [-ceil(thd_lim/4)*4,ceil(thd_lim/4)*4];
[~,ytk.thd] = get_axis_tick(ylm.thd,ylm.thd,ceil(ylm.thd(2)/2),ceil(ylm.thd(2)/4));
ylm.thd = ylm.thd;
ytk.thd = ytk.thd;

ylm.psa = [0,ceil(psa_lim/4)*4];
[~,ytk.psa] = get_axis_tick(ylm.psa,ylm.psa,ceil(ylm.psa(2)/4),ceil(ylm.psa(2)/4));
ylm.psa = ylm.psa;
ytk.psa = ytk.psa;
