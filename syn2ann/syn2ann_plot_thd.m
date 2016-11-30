%% *DISPLACEMENT TIME HISTORY*
max_val = 0;

for i_ = 1:spg(1)
    xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    ypl{i_} = xplot{i_}.syn{identity}.thd.(cpp)*utd.thd;
    max_val = max([max_val;max(abs(ypl{i_}))]);
    leg(i_) = {legplot(i_)};
end
xlbt = cell(spg(1),1);
xlbt(end)= xlb.thd;
ylm.thd = [-ceil(max_val/40)*40,ceil(max_val/40)*40];
[~,ytk.thd] = get_axis_tick(ylm.thd,ylm.thd,ceil(max_val/40)*10,ceil(max_val/40)*10);
ylm.thd = {ylm.thd};
ytk.thd = {ytk.thd};

%
% # _ASIDE_
%
if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.thd.c,'spg',spg,'pax',pax,...
        'scl',scl.thd,'grd',grd.thd,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',xlm.thd(end),'xlb',xlb.thd,'xtk',xtk.thd(end),...
        'ylm',ylm.thd,'ylb',ylb.thd,'ytk',ytk.thd,...
        'leg',leg,'tit',std);
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.thd.c,'spg',spg,'pax',pax,...
        'scl',scl.thd,'grd',grd.thd,...
        'xlm',xlm.thd(end),'xlb',xlb.thd,'xtk',xtk.thd(end),...
        'ylm',ylm.thd,'ylb',ylb.thd,'ytk',ytk.thd,...
        'leg',leg,'tit',std);
end
saveas(gcf,strcat(fnm,sprintf('_thd_%s',cpp)),'epsc');