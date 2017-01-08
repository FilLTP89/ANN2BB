%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _train_ann_justPSA_: function train ANN on PSA values
%% *N.B.*
% Need for:
% _trann_define_inout.m, trann_check_vTn.m,ANN MATLAB tool_
%% *REFERENCES*
% https://fr.mathworks.com/help/nnet/ug/improve-neural-network-generalization-and-avoid-overfitting.html
% http://www.cs.cmu.edu/afs/cs/Web/Groups/AI/html/faqs/ai/neural/faq.html
function train_ann_justPSA(varargin)
    %% *SET-UP*
    wd  = varargin{1};
    ann = varargin{2};
    %
    % _load database_
    %
    db = load(ann.dbn);
    db.nr  = size(db.SIMBAD,2);
    db.vTn = (0:0.05:10)';
    db.nT  = numel(db.vTn);
    %
    % _define input/target natural periods_
    %
    [inp.vTn,tar.vTn,inp.nT,tar.nT] = trann_define_inout(ann.TnC);
    %
    % _check input/target natural periods with database_
    %
    [inp.idx,tar.idx] = trann_check_vTn(inp,tar,db,1e-8);
    %
    % _select class-compatible sites (EC8)_
    %
    switch upper(ann.scl)
        case 'ALL'
            idx_cl = ones(db.nr,1);
        case 'AB'
            ia1 = strcmpi('A',{db.SIMBAD(:).site_EC8});
            ia2 = strcmpi('A*',{db.SIMBAD(:).site_EC8});
            ia  = logical(ia1+ia2);
            ib1 = strcmpi('B',{db.SIMBAD(:).site_EC8});
            ib2 = strcmpi('B*',{db.SIMBAD(:).site_EC8});
            ib  = logical(ib1+ib2);
            idx_cl = logical(ia+ib);
        case 'CD'
            ia1 = strcmpi('C',{db.SIMBAD(:).site_EC8});
            ia2 = strcmpi('C*',{db.SIMBAD(:).site_EC8});
            ia  = logical(ia1+ia2);
            ib1 = strcmpi('D',{db.SIMBAD(:).site_EC8});
            ib2 = strcmpi('D*',{db.SIMBAD(:).site_EC8});
            ib  = logical(ib1+ib2);
            idx_cl = logical(ia+ib);
    end
    idx_cl1 = find(idx_cl==1);
    db.nr   = numel(idx_cl1);
    for i_=1:db.nr
        db.simbad(i_) = db.SIMBAD(idx_cl1(i_));
    end
    
    %% *DEFINE ANN INPUTS/TARGETS (PSA-T*)*
    PSA = -999*ones(db.nr,db.nT);
    switch ann.cp
        % _HORIZONTAL COMPONENT 1_
        case {'h1'}
            for j_ = 1:db.nr
                PSA(j_,:) = db.simbad(j_).psa_h1(:)';
            end
            % _HORIZONTAL COMPONENT 2_
        case {'h2'}
            for j_ = 1:db.nr
                PSA(j_,:) = db.simbad(j_).psa_h2(:)';
            end
        case 'gh'
            for j_ = 1:db.nr
                PSA(j_,:) = geomean([db.simbad(j_).psa_h1(:)';...
                    db.simbad(j_).psa_h2(:)'],1);
            end
            % _VERTICAL COMPONENT_
        case 'ud'
            for j_ = 1:db.nr
                PSA(j_,:) = db.simbad(j_).psa_v(:)';
            end
            
    end
    
    %% *DEFINE INPUT/TARGET PSA POOL (LOG)*
    inp.simbad  = -999*ones(inp.nT,db.nr);
    tar.simbad  = -999*ones(tar.nT,db.nr);
    for i_=1:inp.nT
        inp.simbad(i_,1:db.nr) = log10(PSA(1:db.nr,inp.idx(i_)))';
    end
    for i_=1:tar.nT
        tar.simbad(i_,1:db.nr) = log10(PSA(1:db.nr,tar.idx(i_)))';
    end
    
    %% *DESIGN BASIC ANN*
    dsg = train_ann_basics(ann,db.nr);
    NNs = cell(dsg.ntr,1);
    prf.vld = -999*ones(dsg.ntr,1);
    for i_=1:dsg.ntr
        
        fprintf('ANN %u/%u: \n',i_,dsg.ntr);
        
        %% *DEFINE INPUTS/TARGETS*
        % _ALL INPUT/TARGET TRAINING VALUES_
        NNs{i_}.inp.trn = inp.simbad(:,dsg.idx.trn);
        NNs{i_}.tar.trn = tar.simbad(:,dsg.idx.trn);
        % _ALL INPUT/TARGET VALIDATION VALUES_
        NNs{i_}.inp.vld = inp.simbad(:,dsg.idx.vld);
        NNs{i_}.tar.vld = tar.simbad(:,dsg.idx.vld);
        
        %% *TRAINING ANN*
        % getting net and infos on training sets and performances
        fprintf('TRAINING...\n');
        [NNs{i_}.net,NNs{i_}.trs] = train(dsg.net,NNs{i_}.inp.trn,NNs{i_}.tar.trn);
        
        %% *TEST/VALIDATE ANN PERFORMANCE*
        NNs{i_} = train_ann_valid(NNs{i_});
        prf.vld(i_) = NNs{i_}.prf.vld;
    end
    
    [bst.prf,bst.idx] = min(prf.vld);
    
    %% *PLOT PERFORMANCES*
    % _COMPARE TRAINING_
    set(0,'defaultaxescolororder',[0,0,0;1,0,0;0,0,1;0,1,0]);
    fpplot('xpl',{(1:dsg.ntr)';[1;dsg.ntr];[bst.idx;bst.idx];bst.idx},...
        'ypl',{prf.vld;[mean(prf.vld);mean(prf.vld)];[0;1.01*max(prf.vld)];bst.prf},...
        'pfg',[0,0,15,10],'tit',{'ANN-performance'},...
        'xlb',{'ANN-id'},'xlm',{[0;dsg.ntr+1]},'xtk',{[1;(5:5:dsg.ntr)']},...
        'ylb',{'MSE [1]'},'ylm',{[0;1.05*max(prf.vld)]},'ytk',{[0;1.0*max(prf.vld)]},...
        'lst',{'none';'--';'--';'none'},'lwd',[0.1;2;2;.1],'mrk',{'o';'none';'none';'p'});
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_all_prf_',num2str(dsg.nhn),'n')),'epsc');
    close(gcf);
    
    % _BEST ANN PERFORMANCE_
    set(0,'defaultaxescolororder',[0,0,0]);
    xlm = [NNs{bst.idx}.trs.epoch(1);NNs{bst.idx}.trs.epoch(end)];
    ylm = [1e-3;1e1];
    fpplot(...
        'xpl',{NNs{bst.idx}.trs.epoch';NNs{bst.idx}.trs.epoch';NNs{bst.idx}.trs.epoch';...
        [NNs{bst.idx}.trs.best_epoch;NNs{bst.idx}.trs.best_epoch];xlm},...
        'ypl',{NNs{bst.idx}.trs.perf;NNs{bst.idx}.trs.vperf;NNs{bst.idx}.trs.tperf;
        ylm;...
        [NNs{bst.idx}.trs.vperf(NNs{bst.idx}.trs.best_epoch);...
        NNs{bst.idx}.trs.vperf(NNs{bst.idx}.trs.best_epoch)]},...
        'pfg',[0,0,10,10],'tit',{'ANN-performance'},'scl',{'sly'},...
        'xlb',{'Epochs'},'xlm',{xlm},'leg',{{'TRAIN';'VALID';'TEST'}},...
        'ylb',{'MSE [1]'},'ylm',{ylm},'ytk',{10.^(-3:1)'},...
        'mrk',{'o';'d';'v';'none';'none'},...
        'lst',{'-';'--';':';'-';'-'},'lwd',[2;2;2;1;1]);
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_prf_',num2str(dsg.nhn),'n')),'epsc');
    close(gcf);
    
    % _PLOT REGRESSION_
    % _train set_
    % maximum-minimum values
    col             = jet(size(NNs{bst.idx}.out_trn.trn,1));
    cols            = [col;col];
    set(0,'defaultaxescolororder',cols);
    
    clear xlm ylm
    
    % train set
    xlm.trn = [min(NNs{bst.idx}.out_trn.trn(:));...
        max(NNs{bst.idx}.out_trn.trn(:))];
    ylm.trn = [min(NNs{bst.idx}.out_tar.trn(:));...
        max(NNs{bst.idx}.out_tar.trn(:))];
    xpl.trn = cell(2*size(NNs{bst.idx}.out_trn.trn,1),1);
    ypl.trn = cell(2*size(NNs{bst.idx}.out_trn.trn,1),1);
    % validation set
    xlm.vld = [min(NNs{bst.idx}.out_trn.vld(:));...
        max(NNs{bst.idx}.out_trn.vld(:))];
    ylm.vld = [min(NNs{bst.idx}.out_tar.vld(:));...
        max(NNs{bst.idx}.out_tar.vld(:))];
    xpl.vld = cell(2*size(NNs{bst.idx}.out_trn.vld,1),1);
    ypl.vld = cell(2*size(NNs{bst.idx}.out_trn.vld,1),1);
    % test set
    xlm.tst = [min(NNs{bst.idx}.out_trn.tst(:));...
        max(NNs{bst.idx}.out_trn.tst(:))];
    ylm.tst = [min(NNs{bst.idx}.out_tar.tst(:));...
        max(NNs{bst.idx}.out_tar.tst(:))];
    xpl.tst = cell(2*size(NNs{bst.idx}.out_trn.tst,1),1);
    ypl.tst = cell(2*size(NNs{bst.idx}.out_trn.tst,1),1);
    
    mrk = cell(2*size(NNs{bst.idx}.out_trn.trn,1),1);
    lst = cell(2*size(NNs{bst.idx}.out_trn.trn,1),1);
    lwd = -999*ones(2*size(NNs{bst.idx}.out_trn.trn,1),1);
    leg = cell(2*size(NNs{bst.idx}.out_trn.trn,1),1);
    
    % compute misfit at each period
    for i_=1:size(NNs{bst.idx}.out_trn.trn,1)
        % _TRAIN SET_
        xpl.trn{i_,1} = NNs{bst.idx}.out_trn.trn(i_,:);   % data
        xpl.trn{size(NNs{bst.idx}.out_trn.trn,1)+i_,1} = xlm.trn;                              % linear fit
        ypl.trn{i_,1} = NNs{bst.idx}.out_tar.trn(i_,:);   % data
        % polyfit coefficients
        cfc.trn = polyfit(xpl.trn{i_,1},ypl.trn{i_,1},1);
        % compute linear trend
        ypl.trn{size(NNs{bst.idx}.out_trn.trn,1)+i_,1} = polyval(cfc.trn,xlm.trn);
        % compute Rsquared
        [r2.trn(i_),rmse.trn(i_)] = ...
            rsquare(ypl.trn{i_,1},polyval(cfc.trn,xpl.trn{i_,1}));
        
        % _VALIDATION SET_
        xpl.vld{i_,1} = NNs{bst.idx}.out_trn.vld(i_,:);   % data
        xpl.vld{size(NNs{bst.idx}.out_trn.vld,1)+i_,1} = xlm.vld;                              % linear fit
        ypl.vld{i_,1} = NNs{bst.idx}.out_tar.vld(i_,:);   % data
        % polyfit coefficients
        cfc.vld = polyfit(xpl.vld{i_,1},ypl.vld{i_,1},1);
        % compute linear trend
        ypl.vld{size(NNs{bst.idx}.out_trn.vld,1)+i_,1} = polyval(cfc.vld,xlm.vld);
        % compute Rsquared
        [r2.vld(i_),rmse.vld(i_)] = ...
            rsquare(ypl.vld{i_,1},polyval(cfc.vld,xpl.vld{i_,1}));
        
        % _TEST SET_
        xpl.tst{i_,1} = NNs{bst.idx}.out_trn.tst(i_,:);   % data
        xpl.tst{size(NNs{bst.idx}.out_trn.tst,1)+i_,1} = xlm.tst;                              % linear fit
        ypl.tst{i_,1} = NNs{bst.idx}.out_tar.tst(i_,:);   % data
        % polyfit coefficients
        cfc.tst = polyfit(xpl.tst{i_,1},ypl.tst{i_,1},1);
        % compute linear trend
        ypl.tst{size(NNs{bst.idx}.out_trn.tst,1)+i_,1} = polyval(cfc.tst,xlm.tst);
        % compute Rsquared
        [r2.tst(i_),rmse.tst(i_)] = ...
            rsquare(ypl.tst{i_,1},polyval(cfc.tst,xpl.tst{i_,1}));
        
        mrk{i_,1} = 'o';
        mrk{size(NNs{bst.idx}.out_trn.trn,1)+i_,1} = 'none';
        lst{i_,1} = 'none';
        lst{size(NNs{bst.idx}.out_trn.trn,1)+i_,1} = '-';
        lwd(i_)   = 0.1;
        lwd(size(NNs{bst.idx}.out_trn.trn,1)+i_,1)   = 3.0;
        leg{i_,1} = '';
        leg{size(NNs{bst.idx}.out_trn.trn,1)+i_,1} = ...
            strcat('$T = ',num2str(tar.vTn(i_),'%.2f'),'$');
    end
    xlm.all(1) = floor(min([xlm.trn(1),xlm.vld(1),xlm.tst(1),...
        ylm.trn(1),ylm.vld(1),ylm.tst(1)]));
    xlm.all(2) =  ceil(max([xlm.trn(1),xlm.vld(1),xlm.tst(1),...
        ylm.trn(2),ylm.vld(2),ylm.tst(2)]));
    ylm.all = xlm.all;
    [xtk,ytk] = get_axis_tick(xlm.all,ylm.all,...
        abs(diff(xlm.all))/4,abs(diff(ylm.all))/4);
    
    % _TRAIN SET_
    fpplot('xpl',xpl.trn,'ypl',ypl.trn,'pfg',[0,0,10,10],...
        'mrk',mrk,'lwd',lwd,'lst',lst,'tit',{'Training Set'},...
        'xlb',{'log_10(PSA-ANN)'},'xlm',{xlm.all},'xtk',{xtk},...
        'ylb',{'log_10(PSA-TAR)'},'ylm',{ylm.all},'ytk',{ytk},'leg',{leg});
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_rgr_trn_',num2str(dsg.nhn),'n')),'epsc');
    
    % _VALIDATION SET_
    fpplot('xpl',xpl.vld,'ypl',ypl.vld,'pfg',[0,0,10,10],...
        'mrk',mrk,'lwd',lwd,'lst',lst,'tit',{'Validation Set'},...
        'xlb',{'log_10(PSA-ANN)'},'xlm',{xlm.all},'xtk',{xtk},...
        'ylb',{'log_10(PSA-TAR)'},'ylm',{ylm.all},'ytk',{ytk},'leg',{leg});
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_rgr_vld_',num2str(dsg.nhn),'n')),'epsc');
    
    % _TEST SET_
    fpplot('xpl',xpl.tst,'ypl',ypl.tst,'pfg',[0,0,10,10],...
        'mrk',mrk,'lwd',lwd,'lst',lst,'tit',{'Testing Set'},...
        'xlb',{'log_10(PSA-ANN)'},'xlm',{xlm.all},'xtk',{xtk},...
        'ylb',{'log_10(PSA-TAR)'},'ylm',{ylm.all},'ytk',{ytk},'leg',{leg});
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_rgr_tst_',num2str(dsg.nhn),'n')),'epsc');
    
    set(0,'defaultaxescolororder',[0,0,0]);
    fpplot('xpl',{tar.vTn;tar.vTn;tar.vTn},'ypl',{r2.trn;r2.vld;r2.tst},...
        'pfg',[0,0,10,10],'scl',{'sly'},...
        'mrk',{'o';'d';'v'},'lwd',[1;1;1],'lst',{'-';'--';':'},...
        'xlb',{'T [s]'},'xlm',{[0;tar.vTn(end)]},'xtk',{0:.1:tar.vTn(end)},...
        'ylb',{'R^2 [1]'},'ylm',{[1e-1;1]},'ytk',{(1:10)./10});
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_r2_trn_',num2str(dsg.nhn),'n')),'epsc');
    close all;
    
    
    %% *OUTPUT*
    clear net
    net = NNs{bst.idx};
    save(fullfile(wd,sprintf('net_%u_%s_%s_%un.mat',...
        round(ann.TnC*100),ann.scl,ann.cp,dsg.nhn)),'net');
    return
end