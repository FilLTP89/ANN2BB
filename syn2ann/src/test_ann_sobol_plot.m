plot_set_up;
close all;
%hfg=figure('position',[0,0,12,12]);
xpl = cell(round(srt),1);
ypl = cell(round(srt),1);
leg = cell(round(srt),1);
col =  hsv(round(srt));
%col = [col(1:2:end,:);col(2:2:end,:)];
set(0,'defaultaxescolororder',col);
for k_=1:round(srt)
    xpl{k_,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
    ypl{k_,1} = S{k_,1}(:,2)/2;
    %xpl{2*k_-0,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
    %ypl{2*k_-0,1} = S{k_,1}(:,3)/2.;
    ypl{k_,1}(ypl{k_,1}<=0.1) = -0.0001;
    leg{k_,1} = sprintf('$Sa[Tn = %.1f]$',hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).tid(k_)));
end
[hfg,hax,hpl]=fpplot('xpl',xpl,'ypl',ypl,'pfg',[0 0 12 12],'scl',{'lin'},...
        'xlb',{'T [s]'},'xlm',{[xpl{end,1}(1),xpl{end,1}(end)]},'xtk',{[xpl{end,1}(1):.5:xpl{end,1}(end)]},...
        'ylb',{'[1]'  },'ylm',{[-.01,1.0]                 },'ytk',{0:.25:1},...
        'lwd',1,'lst',{'-'},'mrk',{'o'},'leg',{leg},'tit',{'Sobol Indices-ANN2BB'});
set(hax,'TickLabelInterpreter', 'latex');
saveas(hfg,sprintf('/tmp1/gattif/ann_sobol/results_sobol_ann2bb_%u_%u.eps',i_,j_),'epsc');
close all;


%xpl = cell(srt-round(srt/2),1);
%ypl = cell(srt-round(srt/2),1);
%leg = cell(srt-round(srt/2),1);
%col = hsv(srt-round(srt/2));
%%col = [col(1:2:end,:);col(2:2:end,:)];
%
%for kk_=round(srt/2)+1:srt
%    k_ = kk_-round(srt/2)
%    xpl{k_,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
%    ypl{k_,1} = S{kk_,1}(:,2)/2.;
%    %xpl{2*k_-0,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
%    %ypl{2*k_-0,1} = S{kk_,1}(:,3)/2.;
%    ypl{k_,1}(ypl{k_,1}<=0.1) = -0.0001;
%    leg{k_,1} = sprintf('$Sa[Tn = %.1f]$',hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).tid(kk_)));
%    leg{k_,1} = '';
%end
%[hfg,hax,hpl]=fpplot('xpl',xpl,'ypl',ypl,'pfg',[0 0 12 12],'scl',{'lin'},...
%        'xlb',{'T [s]'},'xlm',{[xpl{end,1}(1),xpl{end,1}(end)]},'xtk',{[xpl{end,1}(1):.5:xpl{end,1}(end)]},...
%        'ylb',{'[1]'  },'ylm',{[-.01,1.0]                 },'ytk',{0:.25:1},...
%        'lwd',2,'lst',{'-'},'leg',{leg},'tit',{'Sobol Indices-ANN2BB'});
%set(hax,'TickLabelInterpreter', 'latex');
%saveas(hfg,sprintf('/tmp1/gattif/ann_sobol/results_sobol_ann2bb_%u_%u_part2.eps',i_,j_),'epsc');
%close all;
