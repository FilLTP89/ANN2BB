%% *VELOCITY TIME HISTORY*
for i_ = 1:spg(1)
    xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    ypl{i_} = xplot{i_}.syn{identity}.thv.(cpp)*utd.thv;
    leg(i_) = {legplot(i_)};
end
%
% # _ASIDE_
%
if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.thv,'grd',grd.thv,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',{xlm.thv},'xlb',xlb.thv,'xtk',{xtk.thv},...
        'ylm',{ylm.thv},'ylb',ylb.thv,'ytk',{ytk.thv},...
        'leg',leg,'tit',std);
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.thv,'grd',grd.thv,...
        'xlm',{xlm.thv},'xlb',xlb.thv,'xtk',{xtk.thv},...
        'ylm',{ylm.thv},'ylb',ylb.thv,'ytk',{ytk.thv},...
        'leg',leg,'tit',std);
end
saveas(gcf,strcat(fnm,sprintf('_thv_%s',cpp)),'epsc');