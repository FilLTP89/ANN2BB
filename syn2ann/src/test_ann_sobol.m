ccc;

global no contr epsilon1
global ann hbs dsx srt
trann_setup_sensitivity_sobol;
syn2ann_pbs_drive;

syn2ann_ann_drive_sobol;

save('VVarEntree.mat','VVarEntree');
dssx = 1:2;
srt = numel(trs.sps.(hbs.mon.cp{1}).tid);
hbs.clc = 2000;
ott = cell(hbs.clc,hbs.mon.na,hbs.mon.nc);
for j_=1:numel(dssx)
    dsx = dssx(j_);
    for i_=1:hbs.mon.na
        load(sprintf('/tmp1/gattif/ann_sobol/results_sobol_ann2bb_%u_%u.mat',i_,j_),'Yy');
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
        for j_=1:hbs.mon.nc
            fnm = sprintf('/tmp1/gattif/ann_sobol/ths_ann2bb_qMC_%s_%s_%u.csv',...
                hbs.mon.cp{j_},hbs.mon.st{i_},k_);
            csvwrite(fnm,spm.sps.(hbs.mon.cp{j_}).syn{i_}.tha.(hbs.mon.cp{j_}));
        end
    end
end

% PLOT SOBOL INDICES
%    for i_=1:hbs.mon.na
%%        [S,Yy,~,~] = Sobol(2,1,0,0,size(VVarEntree{i_,j_},1),hbs.clc,...
%%            VVarEntree{i_,j_},'apply_ann2hbs_sobol(x)',2);
%%        save(sprintf('/tmp1/gattif/ann_sobol/results_sobol_ann2bb_%u_%u.mat',i_,j_),'hbs','trs','S','Yy');
%        load(sprintf('/tmp1/gattif/ann_sobol/results_sobol_ann2bb_%u_%u.mat',i_,j_),'S');
%        test_ann_sobol_plot;
%    end
%end
