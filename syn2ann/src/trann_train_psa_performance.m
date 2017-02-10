function trann_train_psa_performance(varargin)
    %% *SET-UP*
    ann = varargin{1};
    inp = varargin{2};
    tar = varargin{3};
    wd  = varargin{4};
    dsg = varargin{5};
    
    TnC  = inp.vTn(1);
    
    xpl = cell(3,1);
    ypl = cell(3,1);
    
    col = [0,0,0];
    mrk = {'o';'d';'v'};
    lst = {'-';'--';':'};
    leg = {'ALL';'TRAINING';'VALIDATION';'TEST'};
    set(0,'defaultaxescolororder',col);
    
    %% *DEFINE LIMITS*
    % _TRAINING SET_ 
    xlm = [0;1];
    ylm = [-1;1];
    [xtk,ytk] = get_axis_tick(xlm,ylm,abs(diff(xlm))/4,abs(diff(ylm))/4);
    
    
    %% *COMPUTE REGRESSION FOR BEST ANN*
    % _TRAINING SET_ 
    %
    xpl{1,1} = tar.vTn(:)./TnC;
    xpl{2,1} = tar.vTn(:)./TnC;
    xpl{3,1} = tar.vTn(:)./TnC;
    %
    ypl{1,1} = mean(ann.out_trn.all-ann.tar.trn,2);
    err{1,1}(:,1) = abs(ypl{1,1}-min(ann.out_trn.all-ann.tar.trn,[],2));
    err{1,1}(:,2) = abs(ypl{1,1}-max(ann.out_trn.all-ann.tar.trn,[],2));
%     ypl{2,1} = ann.out_trn.vld-ann.out_tar.vld;
%     ypl{3,1} = ann.out_trn.tst-ann.out_tar.tst;
    figure('position',[0,0,10,10]);
    erb = errorbarxy(xpl{1,1},ypl{1,1},...
        zeros(size(xpl{1,1})),zeros(size(xpl{1,1})),...
        err{1,1}(:,1),err{1,1}(:,2),...
        {'ks-', 'k', 'k'});
    erb.hMain.MarkerSize = 15;
    erb.hMain.MarkerEdgeColor = [0,0,0];
    erb.hMain.MarkerFaceColor = [0,0,0];
    keyboard
    errorbar_tick;
    
    xlim(gca,xlm);
    ylim(gca,ylm);
    set(gca,'xtick',xtk,'ytick',ytk);
    xlabel(gca,'T/T_C [1]');
    ylabel(gca,'log_{10}(PSA_{ANN}/PSA_{TAR}) [1]');
    format_figures(gca);
    rule_fig(gcf);
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_lSa_all_',num2str(dsg.nhn),'n')),'epsc');
    
    
    return
end