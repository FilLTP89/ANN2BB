plot_set_up;
close all;
%hfg=figure('position',[0,0,12,12]);
xpl = cell(2*srt,1);
ypl = cell(2*srt,1);
leg = cell(srt,1);
col = parula(srt);

for k_=1:srt
    xpl{2*k_-1,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
    ypl{2*k_-1,1} = S{k_,1}(:,2);
    xpl{2*k_-0,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
    ypl{2*k_-0,1} = S{k_,1}(:,3);
    ypl{2*k_-0,1}(ypl{2*k_-0,1}<=0) = 0.0001;
    leg{k_,1} = sprintf('$Sa[Tn = %.1f]$',hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).tid(k_)));
end
[hfg,hax,hpl]=fpplot('xpl',xpl,'ypl',ypl,'pfg',[0 0 12 12],'scl',{'lin'},...
       'xlb',{'T [s]'},'xlm',{[xpl{1,1}(1),xpl{1,1}(end)]},'xtk',{[xpl{1,1}(1):.5:xpl{1,1}(end)]},...
       'ylb',{'[1]'  },'ylm',{[-.01,1.0]                 },'ytk',{0:.25:1},...
       'lwd',1,'lst',{'--'},'mrk',{'o'},'leg',{leg},'tit',{'Sobol Indices-ANN2BB'});
set(hax,'TickLabelInterpreter', 'latex');
saveas(hfg,sprintf('results_sobol_ann2bb_%u_%u.eps',i_,j_),'epsc');
close all;
