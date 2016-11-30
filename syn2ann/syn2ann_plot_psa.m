%% *SA SPECTRA*
max_val = 0;
for i_ = 1:spg(1)
    xpl{i_} = xplot{i_}.mon.vTn;
    ypl{i_} = abs(xplot{i_}.syn{identity}.psa.(cpp))*utd.psa;
    max_val = max([max_val;max(ypl{i_})]);
end
ylm.psa = [0,ceil(max_val/40)*40];
[~,ytk.psa] = get_axis_tick(ylm.psa,ylm.psa,ceil(max_val/40)*10,ceil(max_val/40)*10);
ylm.psa = {ylm.psa};
ytk.psa = {ytk.psa};
leg = {legplot};

if logical(mod)
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
        'scl',scl.psa,'grd',grd.psa,'mrk',mrkd,'lst',lstd,'lwd',lwdd,...
        'xlm',xlm.psa,'xlb',xlb.psa,'xtk',xtk.psa,...
        'ylm',ylm.psa,'ylb',ylb.psa,'ytk',ytk.psa,...
        'leg',leg,'tit',std);
else
    fpplot('xpl',xpl,'ypl',ypl,'pfg',pfg.psa,...
        'scl',scl.psa,'grd',grd.psa,...
        'xlm',xlm.psa,'xlb',xlb.psa,'xtk',xtk.psa,...
        'ylm',ylm.psa,'ylb',ylb.psa,'ytk',ytk.psa,...
        'leg',leg,'tit',std);
end
saveas(gcf,strcat(fnm,sprintf('_psa_%s',cpp)),'epsc');