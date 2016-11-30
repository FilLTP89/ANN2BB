%% *ACCELERATION TIME HISTORY*
max_val = 0;

for i_ = 1:spg(1)
    xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    ypl{i_} = xplot{i_}.syn{identity}.tha.(cpp)*utd.tha;
    max_val = max([max_val;max(abs(ypl{i_}))]);
    leg(i_) = {legplot(i_)};
end
xlbt = cell(spg(1),1);
xlbt(end)= xlb.tha;
ylm.tha = [-ceil(max_val/40)*40,ceil(max_val/40)*40];
[~,ytk.tha] = get_axis_tick(ylm.tha,ylm.tha,ceil(max_val/40)*10,ceil(max_val/40)*10);
ylm.tha = {ylm.tha};
ytk.tha = {ytk.tha};

%
% # _ASIDE_
%
if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.tha,'grd',grd.tha,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',xlm.tha(end),'xlb',xlb.tha,'xtk',xtk.tha(end),...
        'ylm',ylm.tha,'ylb',ylb.tha,'ytk',ytk.tha,...
        'leg',leg,'tit',std);
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.tha,'grd',grd.tha,...
        'xlm',xlm.tha(end),'xlb',xlb.tha,'xtk',xtk.tha(end),...
        'ylm',ylm.tha,'ylb',ylb.tha,'ytk',ytk.tha,...
        'leg',leg,'tit',std);
end
saveas(gcf,strcat(fnm,sprintf('_tha_%s',cpp)),'epsc');