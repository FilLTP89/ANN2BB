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
psa_lim = 0;
for j_ = 1:numel(cpp.rec)
    psa_lim = max([psa_lim;max(abs(rec.org.syn{mm_}.psa.(cpp.rec{j_})))]);
end
%for j_ = 1:numel(cpp.rec)
%    psa_lim = max([psa_lim;max(abs(hbs.bst.syn{mm_}.psa.(cpp.rec{j_})))]);
%    psa_lim = max([psa_lim;max(abs(spm.sps.(cpp.rec{j_}).syn{mm_}.psa.(cpp.rec{j_})))]);
%end
psa_lim = psa_lim*utd.psa;
ylm.psa = [0,ceil(psa_lim/4)*4];
[~,ytk.psa] = get_axis_tick(ylm.psa,ylm.psa,ceil(ylm.psa(2)/4),ceil(ylm.psa(2)/4));
ylm.psa = ylm.psa;
ytk.psa = ytk.psa;
