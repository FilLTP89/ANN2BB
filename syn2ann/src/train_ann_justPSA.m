%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _train_ann_justPSA_: function train ANN on PSA values
%% *N.B.*
% Need for:
% _trann_define_inout.m, trann_check_vTn.m, trann_tv_sets.m_
% _ANN MATLAB tool_
%% *REFERENCES*
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
    %
    % _define training/validation set_
    %
    [idx_train,idx_valid] = trann_tv_sets(db.nr,5/100);
    
    %% *DEFINE INPUT/OUTPUT*
    PSA = -999*ones(db.nr,db.nT);
    switch ann.cp
        case {'h1'}
            for j_ = 1:db.nr
                PSA(j_,:) = db.simbad(j_).psa_h1(:)';
            end
        case {'h2'}
            for j_ = 1:db.nr
                PSA(j_,:) = db.simbad(j_).psa_h2(:)';
            end
        case 'ud'
            for j_ = 1:db.nr
                PSA(j_,:) = db.simbad(j_).psa_v(:)';
            end
        case 'gh'
            for j_ = 1:db.nr
                PSA(j_,:) = geomean([db.simbad(j_).psa_h1(:)';...
                    db.simbad(j_).psa_h2(:)'],1);
            end
    end
    
    inp.simbad  = -999*ones(inp.nT,db.nr);
    tar.simbad  = -999*ones(tar.nT,db.nr);
    for i_=1:inp.nT
        inp.simbad(i_,1:db.nr) = log10(PSA(1:db.nr,inp.idx(i_)))';
    end
    for i_=1:tar.nT
        tar.simbad(i_,1:db.nr) = log10(PSA(1:db.nr,tar.idx(i_)))';
    end
    
    
    %% *DEFINE NEURON FEATURES*
    % _ACTIVATION FUNCTIONS_
    % func = 'tansig';
    % func = 'purelin';
    % func = 'hardlim';
    % func = 'logsig';
    % The training function BTF can be any of the backprop training
    % functions such as TRAINLM, TRAINBFG, TRAINRP, TRAINGD, etc.
    % *WARNING*: TRAINLM is the default training function because it
    % is very fast, but it requires a lot of memory to run.  If you get
    % an "out-of-memory" error when training try doing one of these:
    %  
    % (1) Slow TRAINLM training, but reduce memory requirements, by
    %     setting NET.efficiency.memoryReduction to 2 or more. (See HELP TRAINLM.)
    % (2) Use TRAINBFG, which is slower but more memory efficient than TRAINLM.
    % (3) Use TRAINRP which is slower but more memory efficient than TRAINBFG.
    %     trainlm is a network training function that updates weight and bias values according to Levenberg-Marquardt optimization.
    %     trainlm is often the fastest backpropagation algorithm in the toolbox, and is highly recommended as a first-choice supervised algorithm, although it does require more memory than other algorithms.
    %     net.trainFcn = 'trainlm' sets the network trainFcn property.
    % Feed-forward networks consist of Nl layers using the DOTPROD
    % weight function, NETSUM net input function, and the specified
    % transfer functions.
    % The first layer has weights coming from the input.  Each subsequent
    % layer has a weight coming from the previous layer.  All layers
    % have biases.  The last layer is the network output.
    % Each layer's weights and biases are initialized with INITNW.
    % Adaption is done with TRAINS which updates weights with the
    % specified learning function. Training is done with the specified
    % training function. Performance is measured according to the specified
    % performance function.
    % There are many types of neural networks.
    % MultiLayerPerceptrons (MLPs) are neural networks with a multiple parallel
    % node-layer topology.
    % Backpropagation (BP) is an algorithm designed for training neural networks
    % with multiple node-layer topologies.
    % There is no such thing as a Backpropagation Network.
    % The obsolete (but still available) functions
    % NEWFIT, NEWPR, NEWFF
    % and the current functions
    % FITNET, PATTERNNET AND FEEDFORWARDNET
    % are all MLPs that, in the default mode, are trained using BP.
    % An alternative time-consuming training approach is to use a
    % genetic algorithm (e.g., GA).
    % NEWFIT (regression and curve-fitting) and NEWPR (classification
    % and pattern-recognition) are specialized algorithms that call
    % NEWFF.
    % FITNET (regression and curve-fitting) and PATTERNNET
    % (classification and pattern-recognition) are specialized algorithms
    % that call FEEDFORWARDNET.
    % The basic difference between FEEDFORWARDNET and FITNET
    % is that the latter yields an additional output: a plot of the output
    % vs target fit.
    % The basic differences between FEEDFORWARDNET and
    % PATTERNNET in the default configurations include
    % 1. Training algorithm: TRAINLM vs TRAINSCG
    % 2. Performance function: MSE vs CROSSENTROPY
    % 3. Output Transfer function: PURELIN vs SOFTMAX
    % 4. Plot Functions
    % All function properties can be obtained by eliminating the ending semicolon
    % in the creation statement. For example:
    % net = patternnet
    
%     % Solve an Input-Output Fitting problem with a Neural Network
%     % Script generated by NFTOOL
%     %
%     % This script assumes these variables are defined:
%     %
%     %   houseInputs - input data.
%     %   houseTargets - target data.
%     
%     inputs = houseInputs;
%     targets = houseTargets;
%     
%     % Create a Fitting Network
%     hiddenLayerSize = 10;
%     net = fitnet(hiddenLayerSize);
%     
%     % Set up Division of Data for Training, Validation, Testing
%     net.divideParam.trainRatio = 70/100;
%     net.divideParam.valRatio = 15/100;
%     net.divideParam.testRatio = 15/100;
%     
%     % Train the Network
%     [net,tr] = train(net,inputs,targets);
%     
%     % Test the Network
%     outputs = net(inputs);
%     errors = gsubtract(outputs,targets);
%     performance = perform(net,targets,outputs)
%     
%     % View the Network
%     view(net)
%     
%     % Plots
%     % Uncomment these lines to enable various plots.
%     % figure, plotperform(tr)
%     % figure, plottrainstate(tr)
%     % figure, plotfit(targets,outputs)
%     % figure, plotregression(targets,outputs)
%     % figure, ploterrhist(errors)
    
    
    %% *CREATE BASE NETWORK - 2LFF-LM (MLP)*
    % Number of Hidden Neurons
    nhn = 30;
    % set up base ANN structure
    % ann.base = newfit(inp.simbad(:,idx_train),tar.simbad(:,idx_train),nhn);
    % ann.base = fitnet(nhn,'trainlm');
    % ann.base = fitnet(nhn,'trainbr');
    ann.base = feedforwardnet(nhn);
    % no show gui
    ann.base.trainParam.showWindow  = false;
    % ann.base.trainParam.show   = 50;
    ann.base.trainParam.lr     = 0.05;
    % Maximum number of epochs to train
    ann.base.trainParam.epochs = 500;
    % performance goal
    ann.base.trainParam.goal   = 1e-3;
    % _percentage of input for training_
    ann.base.divideParam.trainRatio = 85/100;
    % percentage of input for validation
    ann.base.divideParam.valRatio   = 10/100;
    % percentage of input for test
    ann.base.divideParam.testRatio  =  5/100;
    view(ann.base);
    
    %% *TRAIN/TEST ANN NETWORKS*
    % number of trained ANNs
    ann.ntr = 50;
    % ANN-net structure
    ann.net = cell(ann.ntr,1);
    % ANN-training record structure
    ann.trs = cell(ann.ntr,1);
    ann.prf = -999*ones(ann.ntr,1);
    
    for i_=1:ann.ntr
        fprintf('ANN %u/%u: \n',i_,ann.ntr);
        %
        % _TRAINING_
        %
        fprintf('TRAINING...\n');
        [ann.net{i_},ann.trs{i_}] = ...
            train(ann.base,inp.simbad,tar.simbad);
        %
        % _TEST_
        %
        fprintf('TESTING...\n');
        ann.out{i_}.est = ann.net{i_}(inp.simbad);
        ann.out{i_}.trn = ann.out{i_}.est(ann.trs{i_}.trainInd);
        ann.out{i_}.vld = ann.out{i_}.est(ann.trs{i_}.valInd);
        ann.out{i_}.tst = ann.out{i_}.est(ann.trs{i_}.testInd);
        plt.net1 = view(ann.net{i_});
        %
        % _COMPUTE PERFORMANCE_
        %
        fprintf('PERFORMANCE...\n');
        plt.prf = plotperf(ann.trs{i_});
        ann.prf(i_) = perform(ann.net{i_},ann.out{i_},tar.simbad);
        %
        % _REGRESSION_
        %
        ann.tar{i_}.trn = tar.simbad(ann.trs{i_}.trainInd);
        ann.tar{i_}.vld = tar.simbad(ann.trs{i_}.valInd);
        ann.tar{i_}.tst = tar.simbad(ann.trs{i_}.testInd);
    	plt.reg = plotregression(ann.tar{i_}.trn,ann.out{i_}.trn,'Train',...
            ann.tar{i_}.vld,ann.out{i_}.vld,'Validation',...
            ann.tar{i_}.tst,ann.out{i_}.tst,'Testing');
        
        
        
        
        keyboard

    end
    
    % ann.trs{i_} contains all of the information concerning the training 
    % of the network. For example, tr.trainInd, tr.valInd and tr.testInd 
    % contain the indices of the data points that were used in the training, 
    % validation and test sets, respectively. 
    % If you want to retrain the network using the same division of data, 
    % you can set net.divideFcn to 'divideInd', net.divideParam.trainInd to 
    % tr.trainInd, net.divideParam.valInd to tr.valInd, net.divideParam.testInd to tr.testInd.
    % The tr structure also keeps track of several variables during the 
    % course of training, such as the value of the performance function, 
    % the magnitude of the gradient, etc. 
    % You can use the training record to plot the performance progress by 
    % using the plotperf command: plotperf(ann.trs{i_})
    % The property tr.best_epoch indicates the iteration at which the
    % validation performance reached a minimum. The training continued for 
    % 6 more iterations before the training stopped.
    % The validation and test curves are very similar. 
    % If the test curve had increased significantly before the validation 
    % curve increased, then it is possible that some overfitting might have 
    % occurred. The next step in validating the network is to create a 
    % regression plot, which shows the relationship between the outputs of 
    % the network and the targets. If the training were perfect, the network
    % outputs and the targets would be exactly equal, but the relationship
    % is rarely perfect in practice. 
    % For the housing example, we can create a regression plot with the 
    % following commands. 
    % houseOutputs = net(houseInputs);
    % trOut = houseOutputs(tr.trainInd);
    % vOut = houseOutputs(tr.valInd);
    % tsOut = houseOutputs(tr.testInd);
    % trTarg = houseTargets(tr.trainInd);
    % vTarg = houseTargets(tr.valInd);
    % tsTarg = houseTargets(tr.testInd);
    % plotregression(trTarg,trOut,'Train',vTarg,vOut,'Validation',...
    % tsTarg,tsOut,'Testing')
    % The first command calculates the trained network
    % response to all of the inputs in the data set. 
    % The following six commands extract the outputs and targets that 
    % belong to the training, validation and test subsets.
    % The final command creates three regression plots for training, testing and validation.
    
    %% *TEST ANNs AND CHECK BEST PERFORMANCE*
    % Next, each network is tested on the second dataset
    % with both individual performances and the performance for
    % the average output calculated.
    
%     
%     perfs = zeros(1,ann.ntr);
%     output2Total = 0;
    
%     for i_=1:ann.ntr
%         neti         = nets{i_};
%         tri          = tr{i_};
%         output2      = sim(neti,inp.simbad(:,idx_valid));
%         perfs(i_)     = mse(output2-tar.simbad(:,idx_valid));
%         output2Total = output2Total + output2;
% %         plotfit(neti,inp.simbad(:,idx_train),tar.simbad(:,idx_train));
%     end
%     output2AverageOutput = output2Total/ann.ntr;
%     perfAveragedOutputs  = mse(tar.simbad(:,idx_valid)-output2AverageOutput);
%     figure(1)
%     plot([1:ann.ntr],[perfAveragedOutputs].*ones(ann.ntr,1),'--r');
%     % save trained network with the best performance
%     [~,id_min] = min(perfs);
%     net = nets(id_min);
%     net = net{1,1};
    
    %     save(fullfile(wd,sprintf('net_%u_%s_%s_new.mat',...
    %         round(ann.TnC*100),ann.scl,ann.cp)),...
    %         'net','idx_train','idx_valid');
    % Plot
    %         outputs1= sim(net,inp.simbad(:,idx_train));
    
    %         plotregression(targets1(1,:),outputs1(1,:),'1',targets1(2,:),outputs1(2,:),'2',...
    %             targets1(3,:),outputs1(3,:),'3',targets1(4,:),outputs1(4,:),'4',...
    %             targets1(5,:),outputs1(5,:),'5',targets1(6,:),outputs1(6,:),'6');
    return
end

%     vTn_simbad = [0.000000;
%         0.040000;
%         0.040800;
%         0.041700;
%         0.042600;
%         0.043500;
%         0.044400;
%         0.045500;
%         0.046500;
%         0.047600;
%         0.048800;
%         0.050000;
%         0.051300;
%         0.052600;
%         0.054100;
%         0.055600;
%         0.057100;
%         0.058800;
%         0.060600;
%         0.062500;
%         0.064500;
%         0.066700;
%         0.067800;
%         0.069000;
%         0.070200;
%         0.071400;
%         0.072700;
%         0.074100;
%         0.075500;
%         0.076900;
%         0.078400;
%         0.080000;
%         0.081600;
%         0.083300;
%         0.085100;
%         0.087000;
%         0.088900;
%         0.090900;
%         0.093000;
%         0.095200;
%         0.097600;
%         0.100000;
%         0.102000;
%         0.104000;
%         0.106000;
%         0.109000;
%         0.111000;
%         0.114000;
%         0.116000;
%         0.119000;
%         0.122000;
%         0.125000;
%         0.128000;
%         0.132000;
%         0.135000;
%         0.139000;
%         0.143000;
%         0.147000;
%         0.152000;
%         0.156000;
%         0.161000;
%         0.167000;
%         0.172000;
%         0.179000;
%         0.185000;
%         0.192000;
%         0.200000;
%         0.204000;
%         0.208000;
%         0.213000;
%         0.217000;
%         0.222000;
%         0.227000;
%         0.233000;
%         0.238000;
%         0.244000;
%         0.250000;
%         0.256000;
%         0.263000;
%         0.270000;
%         0.278000;
%         0.286000;
%         0.294000;
%         0.303000;
%         0.312000;
%         0.323000;
%         0.333000;
%         0.345000;
%         0.357000;
%         0.370000;
%         0.385000;
%         0.400000;
%         0.417000;
%         0.435000;
%         0.455000;
%         0.476000;
%         0.500000;
%         0.526000;
%         0.556000;
%         0.588000;
%         0.625000;
%         0.667000;
%         0.714000;
%         0.769000;
%         0.833000;
%         0.909000;
%         1.000000;
%         1.050000;
%         1.110000;
%         1.180000;
%         1.250000;
%         1.330000;
%         1.430000;
%         1.540000;
%         1.670000;
%         1.820000;
%         2.000000;
%         2.200000;
%         2.400000;
%         2.600000;
%         2.800000;
%         3.000000;
%         3.200000;
%         3.400000;
%         3.600000;
%         3.800000;
%         4.000000;
%         4.200000;
%         4.400000;
%         4.600000;
%         4.800000;
%         5.000000;
%         5.200000;
%         5.400000;
%         5.600000;
%         5.800000;
%         6.000000;
%         6.200000;
%         6.400000;
%         6.600000;
%         6.800000;
%         7.000000;
%         7.200000;
%         7.400000;
%         7.600000;
%         7.800000;
%         8.000000;
%         8.200000;
%         8.400000;
%         8.600000;
%         8.800000;
%         9.000000;
%         9.200000;
%         9.400000;
%         9.600000;
%         9.800000;
%         10.000000]';

