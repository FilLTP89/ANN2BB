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
%     xpl{1,1} = tar.vTn(:)./TnC;
    xpl{2,1} = tar.vTn(:)./TnC;
    xpl{3,1} = tar.vTn(:)./TnC;
    xpl{4,1} = tar.vTn(:)./TnC;
    %
%     ypl{1,1} = mean(ann.out_trn.all-ann.tar.trn,2);
    ypl{2,1} = mean(ann.out_trn.trn-ann.out_tar.trn,2);
    ypl{3,1} = mean(ann.out_trn.vld-ann.out_tar.vld,2);
    ypl{4,1} = mean(ann.out_trn.tst-ann.out_tar.tst,2);
    %
%     err{1,1}(:,1) = abs(ypl{1,1}-min(ann.out_trn.all-ann.tar.trn,[],2));
%     err{1,1}(:,2) = abs(ypl{1,1}-max(ann.out_trn.all-ann.tar.trn,[],2));
    
    err{2,1}(:,1) = abs(ypl{2,1}-min(ann.out_trn.trn-ann.out_tar.trn,[],2));
    err{2,1}(:,2) = abs(ypl{2,1}-max(ann.out_trn.trn-ann.out_tar.trn,[],2));
    
    err{3,1}(:,1) = abs(ypl{3,1}-min(ann.out_trn.vld-ann.out_tar.vld,[],2));
    err{3,1}(:,2) = abs(ypl{3,1}-max(ann.out_trn.vld-ann.out_tar.vld,[],2));
    
    err{4,1}(:,1) = abs(ypl{4,1}-min(ann.out_trn.tst-ann.out_tar.tst,[],2));
    err{4,1}(:,2) = abs(ypl{4,1}-max(ann.out_trn.tst-ann.out_tar.tst,[],2));
    
    figure('position',[0,0,15,10]);
    
    %
%     erb.all = errorbarxy(xpl{1,1},ypl{1,1},...
%         zeros(size(xpl{1,1})),zeros(size(xpl{1,1})),...
%         err{1,1}(:,1),err{1,1}(:,2),...
%         {'ks--', 'k', 'k'});
%     erb.all.hMain.MarkerSize = 10;
%     erb.all.hMain.MarkerEdgeColor = [0,0,0];
%     erb.all.hMain.MarkerFaceColor = [0,0,0];
%     %
    erb.trn = errorbarxy(xpl{2,1},ypl{2,1},...
        zeros(size(xpl{2,1})),zeros(size(xpl{2,1})),...
        err{2,1}(:,1),err{2,1}(:,2),...
        {'b-', 'b', 'b'});
    erb.trn.hMain.Color = [0.00,0.00,0.00];
    erb.trn.hMain.Marker='none';
    erb.trn.hMain.LineWidth = 3;
%     erb.trn.hMain.MarkerSize = 15;
%     erb.trn.hMain.MarkerEdgeColor = [0,0,0];
%     erb.trn.hMain.MarkerFaceColor = [0,0,0];
    
    for j_=1:size(erb.trn.hErrorbar,2)
        for i_=1:size(erb.trn.hErrorbar,1)
            set(erb.trn.hErrorbar(i_,j_),'color',[0,0,0],...
                'linewidth',3);
            
        end
    end
   
    %
    erb.vld = errorbarxy(xpl{3,1},ypl{3,1},...
        zeros(size(xpl{3,1})),zeros(size(xpl{3,1})),...
        err{3,1}(:,1),err{3,1}(:,2),...
        {'g-', 'g', 'g'});
    erb.vld.hMain.Color = [0.4,0.4,0.4];
    erb.vld.hMain.LineWidth = 3;
    erb.vld.hMain.Marker='none';
%     erb.vld.hMain.MarkerSize = 11;
%     erb.vld.hMain.MarkerEdgeColor = [0.0,0.0,0.0];
%     erb.vld.hMain.MarkerFaceColor = [0.4,0.4,0.4];
    
    for j_=1:size(erb.vld.hErrorbar,2)
        for i_=1:size(erb.vld.hErrorbar,1)
            set(erb.vld.hErrorbar(i_,j_),'color',[0.4,0.4,0.4],...
                'linewidth',3);
            
        end
    end
    %
    erb.tst = errorbarxy(xpl{4,1},ypl{4,1},...
        zeros(size(xpl{4,1})),zeros(size(xpl{4,1})),...
        err{4,1}(:,1),err{4,1}(:,2),...
        {'r-.', 'r', 'r'});
    erb.tst.hMain.Color = [0,0.85,0];
    erb.trn.hMain.Marker='none';
    erb.tst.hMain.LineWidth = 3;
%     erb.tst.hMain.MarkerSize = 7;
%     erb.tst.hMain.MarkerEdgeColor = [0.0,0.0,0.0];
%     erb.tst.hMain.MarkerFaceColor = [0.7,0.7,0.7];
    
    for j_=1:size(erb.tst.hErrorbar,2)
        for i_=1:size(erb.tst.hErrorbar,1)
            set(erb.tst.hErrorbar(i_,j_),'color',[0.0,0.85,0.0],...
                'linewidth',3);
            
        end
    end
    %
    xlim(gca,xlm);
    ylim(gca,ylm);
    set(gca,'xtick',xtk,'ytick',ytk,'linewidth',2);
    set(gca,'ticklength',[.02,.02]);
    xlabel(gca,'T/T* [1]','fontsize',15,'fontweight','bold');
    ylabel(gca,'log_{10}(Sa_{ANN}/Sa_{TAR}) [1]','fontsize',15,'fontweight','bold');
    leg=legend(gca,{'TRN';'VLD';'TST'});
    
    set(leg,'interpreter','latex','location','northeast',...
        'orientation','horizontal','box','off');
    text(0.6,-1.5,sprintf('$N_{n}^{h}$=%u',ann.net.layerWeights{2,1}.size(2)),'parent',gca,...
        'interpreter','latex','fontsize',18)
    rule_fig(gcf);
    
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_lSa_all_',num2str(dsg.nhn),'n')),'epsc');
    
    
    return
end
