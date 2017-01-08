%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _train_ann_basics_: function to design basics ANN
%% *N.B.*
% Need for:_trann_tv_sets.m,ANN MATLAB tool_

function [varargout] = train_ann_basics(varargin)
    %% *SET-UP*
    ann = varargin{1};
    nbs = varargin{2};
    
    %% *CREATE BASE NETWORK - 2LFF-LM (MLP)*
    % ANN name
    dsg.fnm = sprintf('net_%u_%s_%s_dvl',round(ann.TnC*100),ann.scl,ann.cp);
    % Number of Hidden Neurons
    dsg.nhn = 20;
    % set up base ANN structure
    % dsg.net = newfit(inp.simbad(:,idx_train),tar.simbad(:,idx_train),nhn);
    % dsg.net = fitnet(nhn,'trainlm');
    % dsg.net = fitnet(nhn,'trainbr');
    dsg.net = feedforwardnet(dsg.nhn);
    % no show gui
    dsg.net.trainParam.showWindow  = false;
    % NNs.base.trainParam.show   = 50;
    dsg.net.trainParam.lr     = 0.05;
    % Maximum number of epochs to train
    dsg.net.trainParam.epochs = 500;
    % performance goal
    dsg.net.trainParam.goal   = 1e-3;
    % number of trained ANNs
    dsg.ntr = 1;
    % For early stopping, you must be careful not to use an algorithm that
    % converges too rapidly. If you are using a fast algorithm (like trainlm),
    % set the training parameters so that the convergence is relatively slow.
    % For example, set mu to a relatively large value, such as 1, and set
    % mu_dec and mu_inc to values close to 1, such as 0.8 and 1.5, respectively.
    dsg.net.trainParam.mu = 1.0;
    dsg.net.trainParam.mu_dec = 0.8;
    dsg.net.trainParam.mu_inc = 1.5;
    
    %The training functions trainscg and trainbr usually work well with early stopping.
    
    % Set up Division of Data for Training, Validation, Testing
    %     % _percentage of input for training_
    %     dsg.net.divideParam.trainRatio = 85/100;
    %     % percentage of input for validation
    %     dsg.net.divideParam.valRatio   = 10/100;
    %     % percentage of input for test
    %     dsg.net.divideParam.testRatio  =  5/100;
    %
    %     dsg.net.divideParam.trainRatio = 70/100;
    %     dsg.net.divideParam.valRatio = 15/100;
    %     dsg.net.divideParam.testRatio = 15/100;
    [dsg.idx.trn,dsg.idx.vld] = trann_tv_sets(nbs,10/100);
    
    %% *OUTPUT*
    varargout{1} = dsg;
    return
end