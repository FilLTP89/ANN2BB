ccc;

global no contr epsilon1

trann_setup_sensitivity_sobol;
syn2ann_pbs_drive;

syn2ann_ann_drive_sobol;
%trs.sps = apply_ann2hbs_justPSA(pbs.org,ann);
[S,ST,retourTout] = Sobol(2,1,0,0,4,1000,VVarEntree,'apply_ann2hbs_sobol(x)',2);
