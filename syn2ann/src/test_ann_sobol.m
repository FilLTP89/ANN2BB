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
for j_=1:numel(dssx)
    dsx = dssx(j_);
    for i_=1:hbs.mon.na
        [S,Yy,~,~] = Sobol(2,1,0,0,size(VVarEntree{i_,j_},1),hbs.clc,...
            VVarEntree{i_,j_},'apply_ann2hbs_sobol(x)',2);
        save(sprintf('/tmp1/gattif/ann_sobol/results_sobol_ann2bb_%u_%u.mat',i_,j_),'hbs','trs','S','Yy');
        test_ann_sobol_plot;
    end
end
