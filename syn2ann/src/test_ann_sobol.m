ccc;

global no contr epsilon1
global ann hbs dsx srt
trann_setup_sensitivity_sobol;
syn2ann_pbs_drive;

syn2ann_ann_drive_sobol;

save('/tmp1/gattif/ann_sobol_2020/VVarEntree_2020.mat','VVarEntree');
dssx = 1:3;
srt = numel(trs.sps.(hbs.mon.cp{1}).tid);
hbs.clc = 2000; 

flag = 'pp-uq-rho'; %'pp'
flag = 'run';
flag = 'pp-uq-sobol'; 
%flag = 'run';
if strcmpi(flag,'pp-uq-rho')
    s2X = 0.0;
    R   = importdata('/mssmat2/home/gattif/Documents/ares/workdir/ANN2BB/sensitivity/tab8_jea2011.csv'); 
    vTn = R(2:end,1);
    R   = R(2:end,2:end);
    [T,L] = eig(R);
    S = load('/mssmat2/home/gattif/Documents/ares/workdir/ANN2BB/sensitivity/data_rho_fast_u_2019.mat');
    S = S.x.';
    S = uni2norm(S);
    [inp_vTn,tar_vTn,inp_nT,tar_nT] = trann_define_inout(ann.mtd.TnC{1});
    inp.vTn = inp_vTn;
    inp.nT = numel(inp_vTn);
    tar.vTn = tar_vTn;
    tar.nT = numel(tar_vTn);
    [inp_idx,tar_idx] = trann_check_vTn(inp,tar,hbs.mon,1e-8);
    Ri = -999.*ones(numel(inp_vTn),numel(inp_vTn));
    for i_=1:numel(inp_vTn)
        check=find(vTn==inp_vTn(i_));
        if isempty(check) 
            i1 = find(vTn<inp_vTn(i_),1,'last');
            i2 = find(vTn>inp_vTn(i_),1,'first');
            %if abs(vTn(i1)-inp_vTn(i_))>abs(vTn(i2)-inp_vTn(i_))
                tmp = interp1(vTn,R(:,i1),inp_vTn); 
                Ri(i_:end,i_) = [1.;tmp(i_+1:end)];
            %else
            %    tmp = interp1(vTn,R(:,i2),inp_vTn); 
            %    Ri(i_:end,i_) = [1.;tmp(i_:end-1)];
            %end
            % Ri(:,i_) = interp1(vTn,R(:,i1)+(R(:,i2)-R(:,i1))/...
            %     (vTn(i2)-vTn(i1))*(inp_vTn(i_)-vTn(i1)),inp_vTn)
        else
            tmp = interp1(vTn,R(:,check),inp_vTn); 
            Ri(i_:end,i_) = tmp(i_:end); 
        end
        Ri(i_,i_:end) = Ri(i_:end,i_);
    end
    csvwrite('/mssmat2/home/gattif/Documents/ares/workdir/ANN2BB/sensitivity/tab8_jea2011-Ri.csv',...
        [[-1,inp_vTn(:).'];[inp_vTn(:),Ri]]);
    [Ti,Li] = eig(Ri);
    S = Ti*sqrt(Li)*S; 
    hbs.clc = size(S,2); 
    trss = test_ann_sobol_psa2ths_pp_uq_rho(hbs,ann,S,s2X);
    plot_set_up;
    close all;
    %hfg=figure('position',[0,0,12,12]);
    xpl = cell(hbs.mon.nc*(hbs.clc+1),1);
    ypl = cell(hbs.mon.nc*(hbs.clc+1),1);
    leg = cell(hbs.mon.nc*(hbs.clc+1),1);
    col = [1,0,0;0,0,1];
    %col = [col(1:2:end,:);col(2:2:end,:)];
    set(0,'defaultaxescolororder',[repmat(col,[hbs.mon.nc*hbs.clc,1]);0,0,0;0.4,0.4,0.4]);
    for i_ = 1:hbs.mon.na
        for k_=1:hbs.clc
            for j_ = 1:hbs.mon.nc
                xpl{hbs.mon.nc*(k_-1)+j_,1} = [hbs.mon.vTn(trs.sps.(hbs.mon.cp{j_}).tid);
                    hbs.mon.vTn(trs.sps.(hbs.mon.cp{j_}).iid)];
                ypl{hbs.mon.nc*(k_-1)+j_,1} = trss{k_}.(hbs.mon.cp{j_}).syn{i_}.psa.(hbs.mon.cp{j_});
                leg{hbs.mon.nc*(k_-1)+j_,1} = '';
                %xpl{2*k_-0,1} = hbs.mon.vTn(trs.sps.(hbs.mon.cp{dsx}).iid);
                %ypl{2*k_-0,1} = S{k_,1}(:,3)/2.;
            end
        end
        xpl{end-1,1} = hbs.mon.vTn;
        xpl{end+0,1} = hbs.mon.vTn;
        ypl{end-1,1} = hbs.syn{i_}.psa.ew*100;
        ypl{end+0,1} = hbs.syn{i_}.psa.ns*100;
        leg{1,1} = sprintf('%s',upper(hbs.mon.cp{1}));
        leg{2,1} = sprintf('%s',upper(hbs.mon.cp{2}));
        leg{end-1,1} = sprintf('PBS-%s',upper(hbs.mon.cp{1}));
        leg{end-0,1} = sprintf('PBS-%s',upper(hbs.mon.cp{2}));
        [hfg,hax,hpl]=fpplot('xpl',xpl,'ypl',ypl,'pfg',[0 0 12 12],'scl',{'lin'},...
                'xlb',{'T [s]'},'xlm',{[0,3]},'xtk',{0:.5:3},...
                'ylb',{'Sa [cm/s/s]'},'ylm',{[0,400]},'ytk',{0:100:400},...
                'lwd',[ones(hbs.mon.nc*hbs.clc,1);3;3],'lst',{'-'},'leg',{leg},'tit',{sprintf('ANN2BB-%s',hbs.mon.st{i_})});
        set(hax,'TickLabelInterpreter', 'latex');
        saveas(hfg,sprintf('/tmp1/gattif/ann_sobol_2019/sensitivity_ann2bb_%u_%u.eps',i_),'epsc');
        close all;
    end
    %for k_=1:hbs.clc
    %    clc;
    %    disp(sprintf('REALIZATION: %d',k_))
    %    trs = trss{k_};
    %    syn2ann_run;
    %    %% *5). SAVE RESULTS (DNC)*
    %    for i_=1:hbs.mon.na
    %        for j_=1:hbs.mon.nc
    %            fnm = sprintf('/tmp1/gattif/ann_sobol_2019/ths_ann2bb_qMC_%s_%s_%u.csv',...
    %                hbs.mon.cp{j_},hbs.mon.st{i_},k_);
    %            csvwrite(fnm,spm.sps.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}));
    %        end
    %    end
    %end
elseif strcmpi(flag,'pp-uq-sobol')
    hbs.clc = 2000; 
    ott = cell(hbs.clc,hbs.mon.na,hbs.mon.nc);
    for j_=1:numel(dssx)
        dsx = dssx(j_);
        for i_=1:hbs.mon.na
            disp(sprintf('/tmp1/gattif/ann_sobol_2020/results_sobol_ann2bb_%u_%u_2019.mat',i_,j_))
            load(sprintf('/tmp1/gattif/ann_sobol_2020/results_sobol_ann2bb_%u_%u_2019.mat',i_,j_),'Yy');
            for k_=1:hbs.clc
                ott{k_,i_,j_} = cellfun(@(x) x(k_),Yy);
            end
        end
    end
    trss = test_ann_sobol_psa2ths(hbs,ann,ott);

    for k_=1:hbs.clc
        clc;
        disp(sprintf('REALIZATION: %d',k_))
        trs = trss{k_};
        syn2ann_run;
        %% *5). SAVE RESULTS (DNC)*
        for i_=1:hbs.mon.na
            for j_=1:numel(dssx)
                dsx = dssx(j_);
                fnm = sprintf('/tmp1/gattif/ann_sobol_2020/ths_ann2bb_qMC_%s_%s_%u.csv',...
                    hbs.mon.cp{dsx},hbs.mon.st{i_},k_);
                csvwrite(fnm,spm.sps.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}));
            end
        end
    end
elseif strcmpi(flag,'ppsobol')
% PLOT SOBOL INDICES
    for j_=1:numel(dssx)
        dsx = dssx(j_);
        for i_=1:hbs.mon.na
            load(sprintf('/tmp1/gattif/ann_sobol_2019/results_sobol_ann2bb_%u_%u_2019.mat',i_,j_),'S');
            test_ann_sobol_plot;
    end
end
elseif strcmpi(flag,'run')
    for j_=1:numel(dssx)
        dsx = dssx(j_);
        for i_=1:hbs.mon.na
            [S,Yy,~,~] = Sobol(2,1,0,0,size(VVarEntree{i_,j_},1),hbs.clc,...
                VVarEntree{i_,j_},'apply_ann2hbs_sobol(x)',2);
            save(sprintf('/tmp1/gattif/ann_sobol_2020/results_sobol_ann2bb_%u_%u_2019.mat',i_,j_),'hbs','trs','S','Yy');
        end
    end
end
