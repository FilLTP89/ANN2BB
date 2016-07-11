function train_ann(varargin)
    %% *SET-UP*
    ann = varargin{1};
    wd = varargin{2};
    
    database = load(ann.dbn);
    nrec = size(database.SIMBAD,2);
    Tn_simbad = 0:0.05:10;
%     Tn_simbad = [0.000000;
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
    
    nT_simbad = numel(Tn_simbad);
    
    % Define spectral periods for input and output layers of the ANN
    % consider both horizontal components so that nsamples = nrec*2
    % input
    
    if ann.tc == 0.5
        % ANN with corner period = 0.5 s
        Tinput = [0.6:0.1:1.0,1.25:0.25:5.0,-100,100];
        Toutput = [0,0.05,0.1:0.1:0.5];
    elseif ann.tc == 0.6
        % ANN with corner period = 0.6 s
        Tinput = [0.7:0.1:1.0,1.25:0.25:5.0,-100,100];
        Toutput = [0,0.05,0.1:0.1:0.6];
        
    elseif ann.tc == 0.75
        % ANN with corner period = 0.75 s
        Tinput = [0.8:0.1:1.0,1.25:0.25:5.0,-100,100]; % last 2 values = PGV and PGD
        Toutput = [0,0.05,0.1:0.1:0.7,0.75];
        
    elseif ann.tc == 1.0
        % ANN with corner period = 1.0 s
        Tinput = [1.25:0.25:5.0,-100,100]; % last 2 values = PGV and PGD
        Toutput = [0,0.05,0.1:0.1:1.0];
    end
    ninput = length(Tinput);
    moutput = length(Toutput);
    
    % Define input and output data for ANN training based on SIMBAD databased
    % and spectral periods as defined above
    PGA = -999*ones(nrec,1);
    PGV = -999*ones(nrec,1);
    PGD = -999*ones(nrec,1);
    PSA = -999*ones(nrec,nT_simbad);
    switch ann.cp
        case {'h1'}
            for j_ = 1:nrec
                PGA(j_,1) = database.SIMBAD(j_).pga(1);
                PGV(j_,1) = database.SIMBAD(j_).pgv(1);
                PGD(j_,1) = database.SIMBAD(j_).pgd(1);
                PSA(j_,:) = database.SIMBAD(j_).psa_h1(:)';
            end
        case {'h2'}
            for j_ = 1:nrec
                PGA(j_,1) = database.SIMBAD(j_).pga(2);
                PGV(j_,1) = database.SIMBAD(j_).pgv(2);
                PGD(j_,1) = database.SIMBAD(j_).pgd(2);
                PSA(j_,:) = database.SIMBAD(j_).psa_h2(:)';
            end
        case 'ud'
            for j_ = 1:nrec
                PGA(j_,1) = database.SIMBAD(j_).pga(3);
                PGV(j_,1) = database.SIMBAD(j_).pgv(3);
                PGD(j_,1) = database.SIMBAD(j_).pgd(3);
                PSA(j_,:) = database.SIMBAD(j_).psa_v(:)';
            end
        case 'gh'
            for j_ = 1:nrec
                PGA(j_,1) = geomean(database.SIMBAD(j_).pga(1:2));
                PGV(j_,1) = geomean(database.SIMBAD(j_).pgv(1:2));
                PGD(j_,1) = geomean(database.SIMBAD(j_).pgd(1:2));
                PSA(j_,:) = geomean([database.SIMBAD(j_).psa_h1(:)';...
                    database.SIMBAD(j_).psa_h2(:)'],1);
            end
    end
    
%     figure
%     hold all
%     for i_=1:1
%         plot(Tn_simbad,PSA(i_,:));
%     end
%     xlim(gca,[0,4]);xlabel(gca,'T [s]');ylabel(gca,'PSA [cm/s/s]');
%     format_figures(gca);rule_fig(gcf);
    
    iTi = -999*ones(ninput-2,1);
    INPUT_SIMBAD = -999*ones(nrec,ninput-2);
    iTo = -999*ones(moutput,1);
    TARGET_SIMBAD = -999*ones(nrec,moutput);
    
    for k=1:ninput-2
        iTi(k) = find(abs(Tn_simbad-Tinput(k))<1e-6);
        INPUT_SIMBAD(1:nrec,k) = PSA(1:nrec,iTi(k));
    end
    
    INPUT_SIMBAD(1:nrec,k+1) = PGV(1:nrec,1);
    INPUT_SIMBAD(1:nrec,k+2) = PGD(1:nrec,1);
    
    for k=1:moutput
        iTo(k) = find(abs(Tn_simbad-Toutput(k))<1e-6);
        TARGET_SIMBAD(1:nrec,k) = PSA(1:nrec,iTo(k));
    end
    
    %=========================================================================%
    %================== ARTIFICIAL NEURAL NETWORK TRAINING ===================%
    %=========================================================================%
    % Network Inputs and Targets
    inputs = log10(INPUT_SIMBAD)';
    targets = log10(TARGET_SIMBAD)';
    
    % define training set
    ind1 = load(ann.ftr);
    % define validation sets
    ind2 = load(ann.fva);
    
    % training set
    inputs1 = inputs(:,ind1);
    targets1 = targets(:,ind1);
    % validation set
    inputs2 = inputs(:,ind2);
    targets2 = targets(:,ind2);
    
    % Create Network
    numHiddenNeurons = 30;  % Adjust as desired
    net = newfit(inputs1,targets1,numHiddenNeurons);
    net.trainParam.show = 50;
    net.trainParam.lr = 0.05;
    net.trainParam.epochs = 500;
    net.trainParam.goal = 1e-3;
    net.divideParam.trainRatio = 85/100;  % Adjust as desired
    net.divideParam.valRatio = 10/100;  % Adjust as desired
    net.divideParam.testRatio = 5/100;  % Adjust as desired
    
    % numNN neural networks are trained and tested
    numNN = 50;
    nets = cell(1,numNN);
    for i=1:numNN
        disp(['Training ' num2str(i) '/' num2str(numNN)])
        nets{i} = train(net,inputs1,targets1);
    end
    
    
    % Next, each network is tested on the second dataset
    % with both individual performances and the performance for
    % the average output calculated.
    
    perfs = zeros(1,numNN);
    output2Total = 0;
    for i=1:numNN
        neti = nets{i};
        output2 = sim(neti,inputs2);
        perfs(i) = mse(output2-targets2);
        output2Total = output2Total + output2;
    end
    output2AverageOutput = output2Total/numNN;
    perfAveragedOutputs = mse(targets2-output2AverageOutput);
    
    
    %     figure(1)
    %     plot(perfs,'ok');
    %     hold on
    %     plot([1:numNN],[perfAveragedOutputs].*ones(numNN,1),'--r');
    
    % save trained network with the best performance
    [minval,id_min] = min(perfs);
    net = nets(id_min);
    net = net{1,1};
    switch ann.tc
        case 0.5
            save(fullfile(wd,sprintf('net_05s_%s.mat',ann.cp)),'net');
        case 0.6
            save(fullfile(wd,sprintf('net_06s_%s',ann.cp)),'net');
        case 0.75
            save(fullfile(wd,sprintf('net_075s_%s',ann.cp)),'net');
        case 1.0
            save(fullfile(wd,sprintf('net_1s_%s',ann.cp)),'net');
    end
    %     % Plot
    %     outputs1= sim(net,inputs1);
    %     % plotperf(tr);
    %     plotfit(net,inputs1,targets1);
    %     plotregression(targets1(1,:),outputs1(1,:),'1',targets1(2,:),outputs1(2,:),'2',...
    %         targets1(3,:),outputs1(3,:),'3',targets1(4,:),outputs1(4,:),'4',...
    %         targets1(5,:),outputs1(5,:),'5',targets1(6,:),outputs1(6,:),'6');
    return
end
