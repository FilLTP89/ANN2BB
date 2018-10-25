ccc;

global no contr epsilon1
global ann hbs dsx srt
trann_setup_sensitivity_sobol;
syn2ann_pbs_drive;

syn2ann_ann_drive_sobol;
%trs.sps = apply_ann2hbs_justPSA(pbs.org,ann);
dsx = 1;
srt = numel(trs.sps.(hbs.mon.cp{j_}).tid);
hbs.clc = 1000;
[S,ST,retourTout] = Sobol(2,1,0,0,size(VVarEntree{1,1},1),hbs.clc,VVarEntree{1,1},'apply_ann2hbs_sobol(x)',2);
