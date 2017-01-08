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
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_all_perf_',num2str(dsg.nhn),'nhn')),'epsc');
    close(gcf);
    
    % _BEST ANN PERFORMANCE_
    set(0,'defaultaxescolororder',[0,0,1;0,1,0;1,0,0;0,0,0;0,0,0]);
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
        'lst',{'-';'-';'-';'--';'--'},'lwd',[2;2;2;2;2]);
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_perf_',num2str(dsg.nhn),'nhn')),'epsc');
    close(gcf);
    
    % _PLOT REGRESSION_
    % _train set_
    plt.reg.vld = plotregression(NNs{bst.idx}.out_tar.trn,NNs{bst.idx}.out_trn.trn,'Train');
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_regr_trn_',num2str(dsg.nhn),'nhn')),'epsc');
    % validation set
    plt.reg.vld = plotregression(NNs{bst.idx}.out_tar.vld,NNs{bst.idx}.out_trn.vld,'Validation');
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_regr_vld_',num2str(dsg.nhn),'nhn')),'epsc');
    % test set
    plt.reg.tst = plotregression(NNs{bst.idx}.out_tar.tst,NNs{bst.idx}.out_trn.tst,'Test');
    saveas(gcf,fullfile(wd,strcat(dsg.fnm,'_bst_regr_tst_',num2str(dsg.nhn),'nhn')),'epsc');
    
    close all;
    
    clear net
    net = NNs{bst.idx};
    save(fullfile(wd,sprintf('net_%u_%s_%s_%unhn.mat',...
        round(ann.TnC*100),ann.scl,ann.cp,dsg.nhn)),'net');
    return
end