%% *TEST ANN SENSITIVITY*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% Copyright 2018_
%% *NOTES*
% _test_ann_sensitivity_: function to test ANN sensitivity
%% *N.B.*
% Need for:

%% *REFERENCES*
% @Book{Book_Haykin_1999_ANN,
%   Title                    = {{Neural Newtworks: a Comprehensive Foundation}},
%   Author                   = {Haykin, S.},
%   Publisher                = {Prentice-Hall International, Inc.},
%   Year                     = {1999},
% 
%   File                     = {Book_Haykin_1999_ANN.pdf:Book_Haykin_1999_ANN.pdf:PDF}
% }

ccc;
%% *TRAIN SET-UP (CUSTOMIZE)*
trann_setup_sensitivity;

%% *TEST TRAINED ANN (DNC)*
trann_load_sensitivity;

trann_test;

trann_plot_sensitivity;

%% *COMPARE CORRELATED VS DECORRELATED*
clear all;
clc;
wkd = '~/Documents';
res0 = load(fullfile(wkd,'ares/workdir/ANN2BB/sensitivity/res_rho_0.mat'));
res1 = load(fullfile(wkd,'ares/workdir/ANN2BB/sensitivity/res_rho_fast.mat'));

xpl=cell(6,1);
ypl=cell(6,1);
lst=cell(6,1);
mrk=cell(6,1);

for i_=1:res0.rec.org.mon.na
    stb(i_,:,1)=res0.ann.tst{1}.syn{i_}.psa.gh(:).'.*100;
end
for i_=1:res1.rec.org.mon.na
    stb(i_,:,2)=res1.ann.tst{1}.syn{i_}.psa.gh(:).'.*100;
end
xpl(1:3,1)=repmat({res0.ann.tst{1}.mon.vTn},[3,1]);
xpl(4:6,1)=repmat({res1.ann.tst{1}.mon.vTn},[3,1]);
ypl{1,1}=geomean(squeeze(stb(:,:,1)),1);
ypl{2,1}=geomean(squeeze(stb(:,:,1)),1)+std(squeeze(stb(:,:,1)),0,1);
ypl{3,1}=geomean(squeeze(stb(:,:,1)),1)-std(squeeze(stb(:,:,1)),0,1);
ypl{4,1}=geomean(squeeze(stb(:,:,2)),1);
ypl{5,1}=geomean(squeeze(stb(:,:,2)),1)+std(squeeze(stb(:,:,2)),0,1);
ypl{6,1}=geomean(squeeze(stb(:,:,2)),1)-std(squeeze(stb(:,:,2)),0,1);
lst={'-','-.','-.','-','-','-'};

cols=[0.55,0.55,0.55;0.55,0.55,0.55;0.55,0.55,0.55;0,0,0;0,0,0;0,0,0];
set(0,'defaultaxescolororder',cols);
[hfg,hax,~]=fpplot('xpl',xpl,'ypl',ypl,'pfg',[0,0,14,12],'scl',{'slx'},'vfg','on',...
    'xlm',{[0.05,3]},'xtk',{[0.05,0.1,0.75,1,2,3]},'ylm',{[0,320]},'ytk',{[0:40:320]},...
    'lst',lst,'lwd',[3,3,3,1,1,1],'xlb',{'\boldmath$T [s]$'},...
    'ylb',{'\boldmath$Sa \left[cm/s^2\right]$'},...
    'leg',{{'$\left[\mu_{GH}\right]_{\rho=0}$','',...
    '$\left[\mu_{GH}\pm\sigma\right]_{\rho=0}$',...
    '$\left[\mu_{GH}\right]_{\rho\neq0}$','',...
    '$\left[\mu_{GH}\pm\sigma\right]_{\rho\neq0}$'}});
vline(hax(1),0.75,'k--');
hax(1).XLabel.Interpreter='latex';
hax(1).YLabel.Interpreter='latex';
%hax.get
%hax(2).Legend.Interpreter='latex';
%hax(2).Legend.Box='off';
%hax(2).Legend.FontSize=16;
%saveas(hfg,'fig_rho_0','epsc');
print('fig_cmp', '-depsc', '-r500');
