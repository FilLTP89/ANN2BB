function [varargout] = trann_train_best_performance(varargin)
    %% *SET-UP*
    NNs = varargin{1};
    prf = varargin{2};
    dsg = varargin{3};
    wd  = varargin{4};
    

    %% *COMPUTE BEST PERFORMANCE*
    [bst.prf,bst.idx] = min(prf.vld);
    keyboard
    ann = NNs{bst.idx};
    prf.avg = mean(prf.vld);
    prf.std = std(prf.vld);
    
    switch ann.train_strategy
        
        case 'classic'
            prf.mse = mse(ann.net,ann.tar.vld,prf.avg);
        case 'bootstrap'
            prf.mse = mse(ann.net,ann.out_tar.vld,prf.avg);
    end

%     %% *COMPARE TRAINING*
%     set(0,'defaultaxescolororder',[0.65,0.65,0.65;0,0,0;0,0,0;0,0,0]);
%     fpplot('xpl',{(1:dsg.ntr)';bst.idx;[bst.idx;bst.idx];[1;dsg.ntr]},...
%         'ypl',{prf.vld;bst.prf;[0.01;0.1];[prf.avg;prf.avg]},...
%         'pfg',[0,0,10,10],'tit',{'ANN-performance'},'scl',{'sly'},...
%         'xlb',{'ANN-ID'},'xlm',{[0;dsg.ntr+1]},'xtk',{[1;(5:5:dsg.ntr)']},...
%         'ylb',{'MSE [1]'},'ylm',{[0.01;0.1]},'ytk',{(1:10)./100},...
%         'lst',{'none';'none';'--';'--'},'lwd',[0.1;0.1;1;2],'mrk',{'o';'x';'none';'none'});
%     saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_all_prf_',num2str(dsg.nhn),'n')),'epsc');
%     close(gcf);
%     
%     % *BEST ANN PERFORMANCE*
%     set(0,'defaultaxescolororder',[0,0,0]);
%     xlm = [ann.trs.epoch(1);ann.trs.epoch(end)];
%     ylm = [1e-2;1e1];
%     fpplot('xpl',{ann.trs.epoch';ann.trs.epoch';ann.trs.epoch';...
%         [ann.trs.best_epoch;ann.trs.best_epoch];xlm},...
%         'ypl',{ann.trs.perf;ann.trs.vperf;ann.trs.tperf;
%         ylm;...
%         [ann.trs.vperf(ann.trs.best_epoch);...
%         ann.trs.vperf(ann.trs.best_epoch)]},...
%         'pfg',[0,0,10,10],'tit',{'ANN-performance'},'scl',{'sly'},...
%         'xlb',{'Epochs'},'xlm',{xlm},'leg',{{'TRAIN';'VALID';'TEST'}},...
%         'ylb',{'MSE [1]'},'ylm',{ylm},'ytk',{10.^(-2:1)'},...
%         'mrk',{'none';'none';'none';'none';'none'},...
%         'lst',{'-';'--';':';'-';'-'},'lwd',[1;1;1;1;1]);
%     saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_prf_',num2str(dsg.nhn),'n')),'epsc');
%     close(gcf);

    %% *OUTPUT*
    varargout{1} = ann;
    return
end
