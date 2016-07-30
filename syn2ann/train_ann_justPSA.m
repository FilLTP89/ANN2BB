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
    switch upper(ann.cl)
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
        
        db.simbad(i_) = db.SIMBAD(idx_cl1)
    end
    db.simbad = ;
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
    
    % Create Network
    numHiddenNeurons = 30;  % Adjust as desired
    net = newfit(inp.simbad(:,idx_train),tar.simbad(:,idx_train),numHiddenNeurons);
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
        nets{i} = train(net,inp.simbad(:,idx_train),tar.simbad(:,idx_train));
    end
    
    
    % Next, each network is tested on the second dataset
    % with both individual performances and the performance for
    % the average output calculated.
    
    perfs = zeros(1,numNN);
    output2Total = 0;
    for i=1:numNN
        neti         = nets{i};
        output2      = sim(neti,inp.simbad(:,idx_valid));
        perfs(i)     = mse(output2-tar.simbad(:,idx_valid));
        output2Total = output2Total + output2;
    end
%     output2AverageOutput = output2Total/numNN;
%     perfAveragedOutputs  = mse(targets2-output2AverageOutput);
    %     figure(1)
    %     plot(perfs,'ok');
    %     hold on
    %     plot([1:numNN],[perfAveragedOutputs].*ones(numNN,1),'--r');
    
    % save trained network with the best performance
    [~,id_min] = min(perfs);
    net = nets(id_min);
    net = net{1,1};
    save(fullfile(wd,sprintf('net_%u_%s_%s.mat',...
        round(ann.TnC*100),ann.cl,ann.cp)),...
        'net','idx_train','idx_valid');
    %     % Plot
    %     outputs1= sim(net,inputs1);
    %     % plotperf(tr);
    %     plotfit(net,inputs1,targets1);
    %     plotregression(targets1(1,:),outputs1(1,:),'1',targets1(2,:),outputs1(2,:),'2',...
    %         targets1(3,:),outputs1(3,:),'3',targets1(4,:),outputs1(4,:),'4',...
    %         targets1(5,:),outputs1(5,:),'5',targets1(6,:),outputs1(6,:),'6');
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

