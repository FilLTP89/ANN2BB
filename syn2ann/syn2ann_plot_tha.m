%% *ACCELERATION TIME HISTORY*
for i_ = 1:spg(1)
    xpl{i_} = xplot{i_}.mon.vtm{identity}-time_shift(i_);
    ypl{i_} = xplot{i_}.syn{identity}.tha.(cpp)*utd.tha;
    leg(i_) = {legplot(i_)};
end
%
% # _ASIDE_
%
if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.tha,'grd',grd.tha,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',{xlm.tha},'xlb',xlb.tha,'xtk',{xtk.tha},...
        'ylm',{ylm.tha},'ylb',ylb.tha,'ytk',{ytk.tha},...
        'leg',leg,'tit',std);
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.tha.c,'spg',spg,'pax',pax,...
        'scl',scl.tha,'grd',grd.tha,...
        'xlm',{xlm.tha},'xlb',xlb.tha,'xtk',{xtk.tha},...
        'ylm',{ylm.tha},'ylb',ylb.tha,'ytk',{ytk.tha},...
        'leg',leg,'tit',std);
end
saveas(gcf,strcat(fnm,sprintf('_tha_%s',cpp)),'epsc');