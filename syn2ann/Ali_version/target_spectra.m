%% *DEFINITION OF THE TARGET SPECTRA*
    start_inp = find(T==inp_Tn(1));
    end_inp = find(T==inp_Tn(end));
    T =[tar_Tn inp_Tn];

    for i=1:3
   
    PSA = num_sim.PSA_orig(:,i);
    PSA_inp = PSA(start_inp:end_inp);
    PSA_inp = PSA_inp.*100; % cm/s2
    inp = log10(PSA_inp);
    if i == 1||2
         outp = sim(net_h,inp);
         PSA_outp = 10.^(outp);
    elseif i == 3
        outp = sim(net_v,inp);
        PSA_outp = 10.^(outp);
    end
    PSA_target(:,i) = [PSA_outp; PSA_inp]; %cm/s2
    PSA_target(:,i) = PSA_target(:,i);
    end
    
    
