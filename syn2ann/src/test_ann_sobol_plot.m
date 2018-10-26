plot_set_up;
%hfg=figure('position',[0,0,12,12]);
xpl = cell(2*srt,1);
ypl = cell(2*srt,1);
leg = cell(srt,1);
col = jet(srt);

for i_=1:srt
    xpl{2*i_-1,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
    ypl{2*i_-1,1} = S{i_,1}(:,2);
    xpl{2*i_-0,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
    ypl{2*i_-0,1} = S{i_,1}(:,3);
    leg{i_,1} = sprintf('$Sa[Tn = %.1f]$',hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).tid(i_)));
end
fpplot('xpl',xpl,'ypl',ypl,'pfg',[0 0 12 12],...
       'xlb',{'T [s]'},'xlm',{[xpl{1,1}(1),xpl{1,1}(end)]},'xtk',xpl(1,1),...
       'ylb',{'[1]'  },'ylm',{[-.05,1.0]                 },'ytk',{0:.25:1},...
       'leg',{leg},'tit',{'Sobol Indices-ANN2BB'});
set(gca,'TickLabelInterpreter', 'latex');
saveas(hfg,'trial_sobol.eps','epsc');:xa

