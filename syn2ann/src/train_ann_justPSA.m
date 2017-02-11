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
    
    dsg = train_ann_basics(ann,db.nr,ann.train_strategy);
    NNs = cell(dsg.ntr,1);
    prf.vld = -999*ones(dsg.ntr,1);
    out.prf = 0.0;
    
    switch ann.train_strategy
        
        case 'classic'
            
            for i_=1:dsg.ntr
                
                fprintf('ANN %u/%u: \n',i_,dsg.ntr);
                %% *DEFINE INPUTS/TARGETS*
                % _ALL INPUT/TARGET TRAINING VALUES_
                NNs{i_}.inp.trn = inp.simbad(:,dsg.idx.trn);
                NNs{i_}.tar.trn = tar.simbad(:,dsg.idx.trn);
                % _ALL INPUT/TARGET VALIDATION VALUES_
                NNs{i_}.inp.vld = inp.simbad(:,dsg.idx.vld);
                NNs{i_}.tar.vld = tar.simbad(:,dsg.idx.vld);
                NNs{i_}.train_strategy = ann.train_strategy;
                
                %% *TRAINING ANN*
                % getting net and infos on training sets and performances
                fprintf('TRAINING...\n');
                [NNs{i_}.net,NNs{i_}.trs] = ...
                    train(dsg.net,NNs{i_}.inp.trn,NNs{i_}.tar.trn);
                
                %% *TEST/VALIDATE ANN PERFORMANCE*
                NNs{i_} = train_ann_valid(NNs{i_});
                prf.vld(i_) = NNs{i_}.prf.vld;
            end
            
        case 'bootstrap'
            
            for i_=1:dsg.ntr
                
                dsg.net.divideParam.trainInd = ...
                    datasample(dsg.trn_idx,numel(dsg.trn_idx));
                
                fprintf('ANN %u/%u: \n',i_,dsg.ntr);
                %% *DEFINE INPUTS/TARGETS*
                % _ALL INPUT/TARGET TRAINING VALUES_
                NNs{i_}.inp.trn = inp.simbad(:,:);
                NNs{i_}.tar.trn = tar.simbad(:,:);
                NNs{i_}.train_strategy = ann.train_strategy;
                
                %% *TRAINING ANN*
                % getting net and infos on training sets and performances
                fprintf('TRAINING...\n');
                [NNs{i_}.net,NNs{i_}.trs] = ...
                    train(dsg.net,NNs{i_}.inp.trn,NNs{i_}.tar.trn);
                
                %% *TEST/VALIDATE ANN PERFORMANCE*
                NNs{i_} = train_ann_valid(NNs{i_});
                prf.vld(i_) = NNs{i_}.prf.vld;
            end
            
    end
    %% *COMPUTE BEST PERFORMANCE*
    NNs = trann_train_best_performance(NNs,prf,dsg,wd);
    
    %     %% *COMPUTE REGRESSIONS*
    %     trann_train_regression(NNs,dsg,tar.vTn,wd);
    
    %% *COMPUTE REGRESSIONS*
    trann_train_psa_performance(NNs,inp,tar,wd,dsg);
    keyboard
    %% *OUTPUT*
    save(fullfile(wd,sprintf('net_%u_%s_%s_%un.mat',...
        round(ann.TnC*100),ann.scl,ann.cp,dsg.nhn)),'NNs');
    return
end
