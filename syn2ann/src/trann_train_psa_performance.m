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
    xlm =  [0.00;1.00];
    ylm =  [-1.00;1.00];
    xtk =  0.00:0.25:1.00;
    ytk = -1.00:0.25:1.00;
    
    xlm = [-0.05;1];
    cols = gray(6);
    
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
    err{2,1}(:,1) = prctile(ann.out_trn.trn-ann.out_tar.trn,16,2);
    err{2,1}(:,2) = prctile(ann.out_trn.trn-ann.out_tar.trn,84,2);
    
    err{3,1}(:,1) = prctile(ann.out_trn.vld-ann.out_tar.vld,16,2);
    err{3,1}(:,2) = prctile(ann.out_trn.vld-ann.out_tar.vld,84,2);
    %
    err{4,1}(:,1) = prctile(ann.out_trn.tst-ann.out_tar.tst,16,2);
    err{4,1}(:,2) = prctile(ann.out_trn.tst-ann.out_tar.tst,84,2);
    
    
    
    
    figure('position',[0,0,15,10]);
    
    pl11=plot(xpl{2,1},ypl{2,1}); hold all;
%     pl11.LineWidth=4;
%     pl11.Color=rgb('lightgrey');
    pl11.LineStyle='none';
    pl21=bar(xpl{2,1},err{2,1}(:,1)); hold all;
    pl3=bar(xpl{2,1},err{2,1}(:,2)); hold all;
    
    pl21.BarWidth=0.9;
    pl3.BarWidth=0.9;
    
    pl21.FaceColor=rgb('lightgrey');
    pl3.FaceColor=rgb('lightgrey');
    
    pl22=plot(xpl{3,1},ypl{3,1}); hold all;
%     pl22.LineWidth=4;
%     pl22.Color=[0.4,0.4,0.4];
    pl22.LineStyle='none';
    pl23=bar(xpl{3,1},err{3,1}(:,1)); hold all;
    pl3=bar(xpl{3,1},err{3,1}(:,2)); hold all;
    pl23.BarWidth=0.5;
    pl3.BarWidth=0.5;
    pl23.FaceColor=[0.4,0.4,0.4];
    pl3.FaceColor=[0.4,0.4,0.4];
    
    pl33=plot(xpl{4,1},ypl{4,1}); hold all;
    pl33.LineStyle='none';
%     pl33.LineWidth=4;
%     pl33.Color=rgb('black');
    pl24=bar(xpl{4,1},err{4,1}(:,1)); hold all;
    pl3=bar(xpl{4,1},err{4,1}(:,2)); hold all;
    pl24.BarWidth=0.2;
    pl3.BarWidth=0.2;
    pl24.FaceColor=rgb('black');
    pl3.FaceColor=rgb('black');
    
    
    
%     pl1.MarkerSize=13;
%     pl1.MarkerFaceColor=rgb('white');
%     pl1.MarkerEdgeColor=rgb('white');
    
%     pl1.MarkerSize=13;
%     pl1.MarkerFaceColor=[.2,.2,.2];
%     pl1.MarkerEdgeColor=[.2,.2,.2];
    
%     pl1.MarkerSize=13;
%     pl1.MarkerFaceColor=rgb('black');
%     pl1.MarkerEdgeColor=rgb('black');
    
    
    xlim(gca,xlm);
    ylim(gca,ylm);
    set(gca,'xtick',xtk,'ytick',ytk,'linewidth',2);
    set(gca,'ticklength',[.02,.02]);
    xlabel(gca,'T/T*','fontsize',15,'fontweight','bold');
    ylabel(gca,'log_{10}(Sa_{ANN}/Sa_{Obs})','fontsize',15,'fontweight','bold');
    leg=legend(gca,[pl21,pl23,pl24],{'TRN';'VLD';'TST'});
    
    set(leg,'interpreter','latex','location','northeast',...
        'orientation','horizontal','box','off');
    
    text(0.7,-0.7,strcat('$T^\star=$',num2str(TnC,'%.2f'),'$s$'),'parent',gca,...
        'interpreter','latex','fontsize',18)
    rule_fig(gcf);
    
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_lSa_all_',num2str(dsg.nhn),'n')),'epsc');
    
    return
end
