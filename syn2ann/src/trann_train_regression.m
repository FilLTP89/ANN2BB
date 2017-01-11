function trann_train_regression(varargin)
    %% *SET-UP*
    ann = varargin{1};
    dsg = varargin{2};
    vTn = varargin{3};
    wd  = varargin{4};

    %% *DEFINE LIMITS*
    % _TRAINING SET_ 
    xlm.trn = [min(ann.out_trn.trn(:));max(ann.out_trn.trn(:))];
    ylm.trn = [min(ann.out_tar.trn(:));max(ann.out_tar.trn(:))];
    % _VALIDATION SET_
    xlm.vld = [min(ann.out_trn.vld(:));max(ann.out_trn.vld(:))];
    ylm.vld = [min(ann.out_tar.vld(:));max(ann.out_tar.vld(:))];
    % _TEST SET_
    xlm.tst = [min(ann.out_trn.tst(:));max(ann.out_trn.tst(:))];
    ylm.tst = [min(ann.out_tar.tst(:));max(ann.out_tar.tst(:))];
    % _ALL_
    xlm.all(1,1) = floor(min([xlm.trn(1),xlm.vld(1),xlm.tst(1),...
        ylm.trn(1),ylm.vld(1),ylm.tst(1)]));
    xlm.all(2,1) =  ceil(max([xlm.trn(1),xlm.vld(1),xlm.tst(1),...
        ylm.trn(2),ylm.vld(2),ylm.tst(2)]));
    ylm.all = xlm.all;
    [xtk,ytk] = get_axis_tick(xlm.all,ylm.all,...
        abs(diff(xlm.all))/4,abs(diff(ylm.all))/4);

    nTn  = size(ann.out_trn.trn,1);
    col  = [0.65,0.65,0.65;rgb('Orange');0,1,0;0,0,0];

    set(0,'defaultaxescolororder',col);

    %% *COMPUTE REGRESSION FOR BEST ANN*
    % _TRAINING SET_ 
    [rgr.trn.r,rgr.trn.m,rgr.trn.b] = regression(ann.out_tar.trn,ann.out_trn.trn);
    [~,rgr.trn.min] = min(rgr.trn.r);
    [~,rgr.trn.max] = max(rgr.trn.r);
    rgr.trn.xpl{1,1} = ann.out_tar.trn(:);
    rgr.trn.xpl{2,1} = xlm.all;
    rgr.trn.xpl{3,1} = xlm.all;
    rgr.trn.xpl{4,1} = xlm.all;
    rgr.trn.ypl{1,1} = ann.out_trn.trn(:);
    rgr.trn.ypl{2,1} = rgr.trn.b(rgr.trn.min)+rgr.trn.m(rgr.trn.min).*xlm.all;  
    rgr.trn.ypl{3,1} = rgr.trn.b(rgr.trn.max)+rgr.trn.m(rgr.trn.max).*xlm.all;  
    rgr.trn.ypl{4,1} = xlm.all;

    fpplot('xpl',rgr.trn.xpl,'ypl',rgr.trn.ypl,'tit',{'Training'},...
        'pfg',[0,0,12,12],...
	'leg',{{'$TRAIN$';...
	strcat('$R^2_{min}=',num2str(rgr.trn.r(rgr.trn.min)),'- T =',num2str(vTn(rgr.trn.min),'%.2f'),'s$');...
	strcat('$R^2_{max}=',num2str(rgr.trn.r(rgr.trn.max)),'- T =',num2str(vTn(rgr.trn.max),'%.2f'),'s$');...
	''}},...
        'mrk',{'o';'none';'none';'none'},'lwd',[0.1;2.5;2.5;1],'lst',{'none';'-';'-';'--'},...
        'xlb',{'log_{10}(PSA-TAR)'},'xlm',{xlm.all},'xtk',{xtk},...
        'ylb',{'log_{10}(PSA-ANN)'},'ylm',{ylm.all},'ytk',{ytk});
    lgg = get(gcf,'children'); 
    set(lgg(1),'location','northwest');
    ggg=get(lgg(2),'children');
    ggg(4).MarkerEdgeColor = [0,0,0];
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_rgrA_trn_',num2str(dsg.nhn),'n')),'epsc');
    close all;
    
    % _VALIDATION SET_ 
    [rgr.vld.r,rgr.vld.m,rgr.vld.b] = regression(ann.out_tar.vld,ann.out_trn.vld);
    [~,rgr.vld.min] = min(rgr.vld.r);
    [~,rgr.vld.max] = max(rgr.vld.r);
    rgr.vld.xpl{1,1} = ann.out_tar.vld(:);
    rgr.vld.xpl{2,1} = xlm.all;
    rgr.vld.xpl{3,1} = xlm.all;
    rgr.vld.xpl{4,1} = xlm.all;
    rgr.vld.ypl{1,1} = ann.out_trn.vld(:);
    rgr.vld.ypl{2,1} = rgr.vld.b(rgr.vld.min)+rgr.vld.m(rgr.vld.min).*xlm.all;  
    rgr.vld.ypl{3,1} = rgr.vld.b(rgr.vld.max)+rgr.vld.m(rgr.vld.max).*xlm.all;  
    rgr.vld.ypl{4,1} = xlm.all;

    fpplot('xpl',rgr.vld.xpl,'ypl',rgr.vld.ypl,'tit',{'Training'},...
        'pfg',[0,0,12,12],...
	'leg',{{'$TRAIN$';...
	strcat('$R^2_{min}=',num2str(rgr.vld.r(rgr.vld.min)),'- T =',num2str(vTn(rgr.vld.min),'%.2f'),'s$');...
	strcat('$R^2_{max}=',num2str(rgr.vld.r(rgr.vld.max)),'- T =',num2str(vTn(rgr.vld.max),'%.2f'),'s$');...
	''}},...
        'mrk',{'o';'none';'none';'none'},'lwd',[0.1;2.5;2.5;1],'lst',{'none';'-';'-';'--'},...
        'xlb',{'log_{10}(PSA-TAR)'},'xlm',{xlm.all},'xtk',{xtk},...
        'ylb',{'log_{10}(PSA-ANN)'},'ylm',{ylm.all},'ytk',{ytk});
    lgg = get(gcf,'children'); 
    set(lgg(1),'location','northwest');
    ggg=get(lgg(2),'children');
    ggg(4).MarkerEdgeColor = [0,0,0];
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_rgrA_vld_',num2str(dsg.nhn),'n')),'epsc');
    close all;

    % _TEST SET_ 
    [rgr.tst.r,rgr.tst.m,rgr.tst.b] = regression(ann.out_tar.tst,ann.out_trn.tst);
    [~,rgr.tst.min] = min(rgr.tst.r);
    [~,rgr.tst.max] = max(rgr.tst.r);
    rgr.tst.xpl{1,1} = ann.out_tar.tst(:);
    rgr.tst.xpl{2,1} = xlm.all;
    rgr.tst.xpl{3,1} = xlm.all;
    rgr.tst.xpl{4,1} = xlm.all;
    rgr.tst.ypl{1,1} = ann.out_trn.tst(:);
    rgr.tst.ypl{2,1} = rgr.tst.b(rgr.tst.min)+rgr.tst.m(rgr.tst.min).*xlm.all;  
    rgr.tst.ypl{3,1} = rgr.tst.b(rgr.tst.max)+rgr.tst.m(rgr.tst.max).*xlm.all;  
    rgr.tst.ypl{4,1} = xlm.all;

    fpplot('xpl',rgr.tst.xpl,'ypl',rgr.tst.ypl,'tit',{'Training'},...
        'pfg',[0,0,12,12],...
	'leg',{{'$TRAIN$';...
	strcat('$R^2_{min}=',num2str(rgr.tst.r(rgr.tst.min)),'- T =',num2str(vTn(rgr.tst.min),'%.2f'),'s$');...
	strcat('$R^2_{max}=',num2str(rgr.tst.r(rgr.tst.max)),'- T =',num2str(vTn(rgr.tst.max),'%.2f'),'s$');...
	''}},...
        'mrk',{'o';'none';'none';'none'},'lwd',[0.1;2.5;2.5;1],'lst',{'none';'-';'-';'--'},...
        'xlb',{'log_{10}(PSA-TAR)'},'xlm',{xlm.all},'xtk',{xtk},...
        'ylb',{'log_{10}(PSA-ANN)'},'ylm',{ylm.all},'ytk',{ytk});
    lgg = get(gcf,'children'); 
    set(lgg(1),'location','northwest');
    ggg=get(lgg(2),'children');
    ggg(4).MarkerEdgeColor = [0,0,0];
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_rgrA_tst_',num2str(dsg.nhn),'n')),'epsc');
    close all;
    
    return
end
