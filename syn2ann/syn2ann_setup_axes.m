pfg.psa = [0 0 12 12];
xlm.psa = {[0;3]};
xtk.psa = {(0:.5:3)'};
xlb.psa = {'T [s]'};
ylb.psa = {'Sa [cm/s/s]','',''};
scl.psa = {'lin','lin','lin'};
grd.psa = {'minor'};
utd.psa = 100;
ttt = strcat(bhr.nm{mm_},{' '},evt2tit(bhr.st{mm_}.ev{1},bhr.st{mm_}.tp{1}));
ttt = ttt{1};

fprintf('PLOTTING %s\n\n',ttt)

switch  ttt
    case 'ACC 2016-10-26 17:10'
        % *======================== ACC ==============================*
        %             vtm_shift(:) = 0;
        vtm_lim = [0;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;320]};
        ytk.psa = {(0:80:320)'};
    case 'ACC 2016-10-26 19:18'
        % *======================== ACC ==============================*
        %             vtm_shift(:) = 0;
        vtm_lim = [0;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;320]};
        ytk.psa = {(0:80:320)'};
    case 'ACC 2016-10-30 06:40'
        % *======================== ACC ==============================*
        %             vtm_shift(:) = 0;
        vtm_lim = [0;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-60;60];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;1600]};
        ytk.psa = {(0:400:1600)'};
    case 'NRC 2016-08-24 01:36'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;2000]};
        ytk.psa = {(0:500:2000)'};
    case 'NRC 2016-08-24 02:33'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;800]};
        ytk.psa = {(0:200:800)'};
    case 'NRC 2016-10-26 17:10'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;800]};
        ytk.psa = {(0:200:800)'};
    case 'NRC 2016-10-26 19:18'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;1600]};
        ytk.psa = {(0:400:1600)'};
    case 'NRC 2016-10-30 06:40'
        % *======================== NRC ============================*
        %             vtm_shift(:) = 0.73;
        vtm_lim = [5;20];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-4e2;4e2];
        %             thv_lim = [-50;50];
        %             thd_lim = [-20;20];
        ylm.psa = {[0;2000]};
        ytk.psa = {(0:500:2000)'};
        
    case 'KMMH16-1 2016-04-14 21:26'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;1200]};
        ytk.psa = {(0:400:1200)'};
    case 'KMMH16-2 2016-04-14 21:26'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3600]};
        ytk.psa = {(0:900:3600)'};
    case 'KMMH16_1 2016-04-14 22:07'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;640]};
        ytk.psa = {(0:160:640)'};
    case 'KMMH16_2 2016-04-14 22:07'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;2000]};
        ytk.psa = {(0:500:2000)'};
    case 'KMMH16_1 2016-04-15 00:03'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;400]};
        ytk.psa = {(0:100:400)'};
    case 'KMMH16_2 2016-04-15 00:03'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;2000]};
        ytk.psa = {(0:500:2000)'};
    case 'KMMH16_1 2016-04-16 01:25'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;800]};
        ytk.psa = {(0:200:800)'};
    case 'KMMH16_2 2016-04-16 01:25'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3600]};
        ytk.psa = {(0:900:3600)'};
    case 'KMMH16_1 2016-04-16 01:46'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;600]};
        ytk.psa = {(0:200:600)'};
    case 'KMMH16_2 2016-04-16 01:46'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3000]};
        ytk.psa = {(0:600:3000)'};
    case 'KMMH14_1 2016-04-14 21:26'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;900]};
        ytk.psa = {(0:300:900)'};
    case 'KMMH14_2 2016-04-14 21:26'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3000]};
        ytk.psa = {(0:600:3000)'};
    case 'KMMH14_1 2016-04-15 00:03'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;600]};
        ytk.psa = {(0:200:600)'};
    case 'KMMH14_2 2016-04-15 00:03'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;3000]};
        ytk.psa = {(0:600:3000)'};
    case 'KMMH14_1 2016-04-16 01:25'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [15;35];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-3e2;3e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;1200]};
        ytk.psa = {(0:400:1200)'};
    case 'KMMH14_2 2016-04-16 01:25'
        % *======================== KMM ==============================*
        %             vtm_shift(:) = 1.25;
        vtm_lim = [12;52];
        vtm_lab = (vtm_lim(1):5:vtm_lim(end));
        tha_lim = [-10e2;10e2];
        %             thv_lim = [-30;30];
        %             thd_lim = [-30;30];
        ylm.psa = {[0;1600]};
        ytk.psa = {(0:400:1600)'};
end