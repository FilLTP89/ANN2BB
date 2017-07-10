function trann_train_psa_performance(varargin)
    %% *SET-UP*
    ann = varargin{1};
    inp = varargin{2};
    tar = varargin{3};
    wd  = varargin{4};
    dsg = varargin{5};
    
    plot_set_up;
    
    TnC  = inp.vTn(1);
    
    xpl = cell(3,1);
    ypl = cell(3,1);
    
    
    %% *DEFINE LIMITS*
    % _TRAINING SET_ 
    xlm = [0;1];
    ylm = [-2;2];
    [xtk,ytk] = get_axis_tick(xlm,ylm,abs(diff(xlm))/4,abs(diff(ylm))/4);
    
    xlm = [-0.03;1];
    %% *COMPUTE ERROR BARS*
    % _TRAINING SET_ 
    %
    xpl{2,1} = tar.vTn(:)./TnC;
    xpl{3,1} = tar.vTn(:)./TnC;
    xpl{4,1} = tar.vTn(:)./TnC;
    %
    ypl{2,1} = mean(ann.out_trn.trn-ann.out_tar.trn,2);
    ypl{3,1} = mean(ann.out_trn.vld-ann.out_tar.vld,2);
    ypl{4,1} = mean(ann.out_trn.tst-ann.out_tar.tst,2);
    %
    err{2,1}(:,1) = abs(ypl{2,1}-min(ann.out_trn.trn-ann.out_tar.trn,[],2));
    err{2,1}(:,2) = abs(ypl{2,1}-max(ann.out_trn.trn-ann.out_tar.trn,[],2));
    
    err{3,1}(:,1) = abs(ypl{3,1}-min(ann.out_trn.vld-ann.out_tar.vld,[],2));
    err{3,1}(:,2) = abs(ypl{3,1}-max(ann.out_trn.vld-ann.out_tar.vld,[],2));
    %
    err{4,1}(:,1) = abs(ypl{4,1}-min(ann.out_trn.tst-ann.out_tar.tst,[],2));
    err{4,1}(:,2) = abs(ypl{4,1}-max(ann.out_trn.tst-ann.out_tar.tst,[],2));
    
    figure('position',[0,0,15,10]);
    erb.trn = errorbarxy(xpl{2,1},ypl{2,1},...
        zeros(size(xpl{2,1})),zeros(size(xpl{2,1})),...
        err{2,1}(:,1),err{2,1}(:,2),...
        {'b-', 'b', 'b'});
    erb.trn.hMain.Color = rgb2gray(rgb('green'));
    erb.trn.hMain.Marker='none';
    erb.trn.hMain.LineWidth = 8;
    
    for j_=1:size(erb.trn.hErrorbar,2)
        for i_=1:size(erb.trn.hErrorbar,1)
            set(erb.trn.hErrorbar(i_,j_),'color',rgb2gray(rgb('green')),...
                'linewidth',8);
        end
    end
   
    %
    erb.vld = errorbarxy(xpl{3,1},ypl{3,1},...
        zeros(size(xpl{3,1})),zeros(size(xpl{3,1})),...
        err{3,1}(:,1),err{3,1}(:,2),...
        {'g-', 'g', 'g'});
    erb.vld.hMain.Color = rgb2gray([0.65,0.65,0.65]);
    erb.vld.hMain.LineWidth = 6;
    erb.vld.hMain.Marker='none';
    
    for j_=1:size(erb.vld.hErrorbar,2)
        for i_=1:size(erb.vld.hErrorbar,1)
            set(erb.vld.hErrorbar(i_,j_),'color',rgb2gray([0.65,0.65,0.65]),...
                'linewidth',6);
            
        end
    end
    %
    erb.tst = errorbarxy(xpl{4,1},ypl{4,1},...
        zeros(size(xpl{4,1})),zeros(size(xpl{4,1})),...
        err{4,1}(:,1),err{4,1}(:,2),...
        {'r-', 'r', 'r'});
    erb.tst.hMain.Color = rgb2gray([0,0,0]);
    erb.trn.hMain.Marker='none';
    erb.tst.hMain.LineWidth = 4;
    
    for j_=1:size(erb.tst.hErrorbar,2)
        for i_=1:size(erb.tst.hErrorbar,1)
            set(erb.tst.hErrorbar(i_,j_),'color',rgb2gray([0,0,0]),...
                'linewidth',4);
            
        end
    end
    %
    xlim(gca,xlm);
    ylim(gca,ylm);
    set(gca,'xtick',xtk,'ytick',ytk,'linewidth',2);
    set(gca,'ticklength',[.02,.02]);
    xlabel(gca,'T/T*','fontsize',15,'fontweight','bold');
    ylabel(gca,'log_{10}(Sa_{ANN}/Sa_{Obs})','fontsize',15,'fontweight','bold');
    leg=legend(gca,{'TRN';'VLD';'TST'});
    
    set(leg,'interpreter','latex','location','northeast',...
        'orientation','horizontal','box','off');
    
    text(0.6,-1.5,strcat('$T^\star$=',num2str(TnC,'%.2f')),'parent',gca,...
        'interpreter','latex','fontsize',18)
    rule_fig(gcf);
    
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_lSa_all_',num2str(dsg.nhn),'n')),'epsc');
    
    return
end
