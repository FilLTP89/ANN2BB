pfg.psa = [0 0 14 14];
xlm.psa = {[0;3]};
xtk.psa = {(0:.5:3)'};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
utd.psa = 100;
switch  bhr.nm{mm_}
    case 'AMT'
        % *======================== AMT ==============================*
        disp('PLOTTING AMT')
        %             vtm_shift(:) = 0;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;750]};
        ytk.psa = {(0:250:750)'};
    case 'NRC'
        % *======================== NRC ============================*
        disp('PLOTTING NRC')
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;1800]};
        ytk.psa = {(0:600:1800)'};
    case {'KMM1'}
        % *======================== KMM ==============================*
        disp('PLOTTING KMMi1')
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;1200]};
        ytk.psa = {(0:400:1200)'};
    case {'KMM2'}
        % *======================== KMM ==============================*
        disp('PLOTTING KMM')
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3000]};
        ytk.psa = {(0:600:3000)'};
end

