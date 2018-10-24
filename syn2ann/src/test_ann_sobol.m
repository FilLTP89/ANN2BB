ccc;

global no contr epsilon1


trann_setup_sensitivity;



[S,ST,retourTout] = Sobol(2,1,0,0,4,1000,VVarEntree,'trann_ann_sensitivity_sobol(x)',2);
