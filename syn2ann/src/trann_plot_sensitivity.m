load(res.fnm);
xpl=cell(2*rec.org.mon.na,1);
ypl=cell(2*rec.org.mon.na,1);
lst=cell(2*rec.org.mon.na,1);
mrk=cell(2*rec.org.mon.na,1);

for i_=1:rec.org.mon.na
    xpl{2*i_-1,1}=ann.tst{1}.mon.vTn;
    xpl{2*i_-1,1}(1)=0.01;
    xpl{2*i_-0,1}=rec.org.mon.vTn(rec.org.idx_sns);
    ypl{2*i_-1,1}=ann.tst{1}.syn{i_}.psa.gh.*100;
    ypl{2*i_-0,1}=rec.org.syn{i_}.psa.ew(rec.org.idx_sns,1).*100;
    lst{2*i_-1,1}='-';
    lst{2*i_-0,1}='-';
    mrk{2*i_-1,1}='none';
    mrk{2*i_-0,1}='none';
    stb(i_,:)=ann.tst{1}.syn{i_}.psa.gh(:).'.*100;
end
cols=repmat([0.6,0.6,0.6],[rec.org.mon.na,1]);
set(0,'defaultaxescolororder',cols);
[hfg,hax,~]=fpplot('xpl',xpl,'ypl',ypl,'pfg',[0,0,14,12],'scl',{'slx'},'vfg','on',...
    'xlm',{[0.05,3]},'xtk',{[0.05,0.1,1,2,3]},'ylm',{[0,320]},'ytk',{[0:40:320]},...
    'lst',lst,'mrk',mrk,'lwd',2,'xlb',{'\boldmath$T [s]$'},...
    'ylb',{'\boldmath$Sa \left[cm/s^2\right]$'});
plgm=plot(ann.tst{1}.mon.vTn,geomean(stb,1),'k','linewidth',3);
plsp=plot(ann.tst{1}.mon.vTn,geomean(stb,1)+std(stb,0,1),'k--','linewidth',3);
plsm=plot(ann.tst{1}.mon.vTn,geomean(stb,1)-std(stb,0,1),'k--','linewidth',3);
hax(1).XLabel.Interpreter='latex';
hax(1).YLabel.Interpreter='latex';

leg=legend(hax(1),[plgm,plsp],{'$\mu_{GH}$','$\mu_{GH}\pm\sigma$'});
leg.Interpreter='latex';
leg.Box='off';
leg.FontSize=16;
%saveas(hfg,'fig_rho_0','epsc');
print('fig_rho_0', '-depsc', '-r500');
