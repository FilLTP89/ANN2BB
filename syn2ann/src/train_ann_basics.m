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
    % _number of Hidden Neurons_
    dsg.nhn = 50;
    % _number of trained ANNs_
    dsg.ntr = 50;
    % _select input records_
    [dsg.idx.trn,dsg.idx.vld] = trann_tv_sets(nbs,10/100);
    % _set up base ANN structure_
    
    dsg.net = feedforwardnet(dsg.nhn,'trainlm');
    % dsg.net = feedforwardnet(dsg.nhn,'trainbfg');
    % dsg.net = feedforwardnet(dsg.nhn,'trainlm');
    % _show/hide gui_
    dsg.net.trainParam.showWindow  = false;
    % _learning rate_
    % dsg.net.trainParam.lr     = 0.05;
    % _maximum number of epochs to train_
    dsg.net.trainParam.epochs = 500;
    % _performance goal_
    dsg.net.trainParam.goal   = 1e-5;
    
    switch dsg.net.trainFcn
        case 'trainlm'
            % dsg.net.trainParam.mu = 1.0;
            % dsg.net.trainParam.mu_dec = 0.8;
            % dsg.net.trainParam.mu_inc = 1.5;
        case 'trainbfg'
            dsg.net.performParam.regularization = 0.5;
    end
    % For early stopping, you must be careful not to use an algorithm that
    % converges too rapidly. If you are using a fast algorithm (like trainlm),
    % set the training parameters so that the convergence is relatively slow.
    % For example, set mu to a relatively large value, such as 1, and set
    % mu_dec and mu_inc to values close to 1, such as 0.8 and 1.5, respectively.
    % The training functions trainscg and trainbr usually work well with early stopping.
    
    % Set up Division of Data for Training, Validation, Testing
    % _percentage of input for training_
    dsg.net.divideParam.trainRatio = 85/100;
    % dsg.net.divideParam.trainRatio = 70/100;
    % _percentage of input for validation_
    dsg.net.divideParam.valRatio   = 10/100;
    % dsg.net.divideParam.valRatio = 15/100;
    % _percentage of input for test_
    dsg.net.divideParam.testRatio  =  5/100;
    % dsg.net.divideParam.testRatio = 15/100;
    
    %% *OUTPUT*
    varargout{1} = dsg;
    return
end
