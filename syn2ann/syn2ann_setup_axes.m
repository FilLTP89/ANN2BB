pfg.psa = [0 0 14 14];
xlm.psa = {[0;3]};
xtk.psa = {(0:.5:3)'};
xlb.psa = {'T [s]'};
ylb.psa = {'PSA [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
utd.psa = 100;
fprintf('PLOTTING %s\n\n',strcat(bhr.nm{mm_},'_',bhr.st{mm_}.ev{1}))
switch  strcat(bhr.nm{mm_},'_',bhr.st{mm_}.ev{1})
    case 'ACC_20161026_171036'
        % *======================== ACC ==============================*
        %             vtm_shift(:) = 0;
        vtm_lim = [0;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;500]};
        ytk.psa = {(0:100:500)'};
    case 'ACC_20161026_191806'
        % *======================== ACC ==============================*
        %             vtm_shift(:) = 0;
        vtm_lim = [0;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;500]};
        ytk.psa = {(0:100:500)'};
    case 'ACC_20161030_064018'
        % *======================== ACC ==============================*
        %             vtm_shift(:) = 0;
        vtm_lim = [0;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;1600]};
        ytk.psa = {(0:400:1600)'};
    case 'NRC_20160824_013632'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;2000]};
        ytk.psa = {(0:500:2000)'};
    case 'NRC_20160824_023329'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;600]};
        ytk.psa = {(0:200:600)'};
    case 'NRC_20161026_171036'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;900]};
        ytk.psa = {(0:300:900)'};
    case 'NRC_20161026_191806'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;1600]};
        ytk.psa = {(0:400:1600)'};
    case 'NRC_20161030_064018'
        % *======================== NRC ============================*
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
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3000]};
        ytk.psa = {(0:600:3000)'};
end

