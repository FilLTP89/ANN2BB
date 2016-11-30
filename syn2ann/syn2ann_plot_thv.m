%% *VELOCITY TIME HISTORY*
max_val = 0;

for i_ = 1:spg(1)
    xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    ypl{i_} = xplot{i_}.syn{identity}.thv.(cpp)*utd.thv;
    max_val = max([max_val;max(abs(ypl{i_}))]);
    leg(i_) = {legplot(i_)};
end
xlbt = cell(spg(1),1);
xlbt(end)= xlb.thv;
ylm.thv = [-ceil(max_val/40)*40,ceil(max_val/40)*40];
[~,ytk.thv] = get_axis_tick(ylm.thv,ylm.thv,ceil(max_val/40)*10,ceil(max_val/40)*10);
ylm.thv = {ylm.thv};
ytk.thv = {ytk.thv};

%
% # _ASIDE_
%
if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.thv,'grd',grd.thv,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',xlm.thv(end),'xlb',xlb.thv,'xtk',xtk.thv(end),...
        'ylm',ylm.thv,'ylb',ylb.thv,'ytk',ytk.thv,...
        'leg',leg,'tit',std);
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.thv,'grd',grd.thv,...
        'xlm',xlm.thv(end),'xlb',xlb.thv,'xtk',xtk.thv(end),...
        'ylm',ylm.thv,'ylb',ylb.thv,'ytk',ytk.thv,...
        'leg',leg,'tit',std);
end
saveas(gcf,strcat(fnm,sprintf('_thv_%s',cpp)),'epsc');