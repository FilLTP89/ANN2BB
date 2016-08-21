%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');
% ============ DEBUG CODE---> CHECK WITH BBSYNT ============================
[BB_DTIME,BB_ACC,BB_FREQ,BB_AMP] = bbsynt_original(nss.org.mon.dtm(1),nss.org.syn{1}.tha.e,...
    sps.org.mon.dtm(1),sps.org.syn{1}.tha.e,mon.fa(1),mon.fb(1));
BB_TIME = BB_DTIME*(0:numel(BB_ACC)-1);
BB_VEL = cumtrapz(BB_ACC)*BB_DTIME;
BB_DIS = cumtrapz(BB_VEL)*BB_DTIME;
BB_FREQ(end+1:numel(BB_AMP)) = BB_FREQ(end)+mean(diff(BB_FREQ))*(1:numel(BB_AMP)/2)';

%% *RESAMPLING*
fprintf('--> Resampling\n');
[nss.org,sps.org] = lfhf_rsmpl(nss.org,sps.org);

%% *PAD LF/HF*
fprintf('--> Padding/Tapering\n');
[nss.org,sps.org] = lfhf_pad(nss.org,sps.org);

%% *ALIGN LF/HF*
fprintf('--> Align records\n');
[nss.org,sps.org] = lfhf_shift(nss.org,sps.org);
nss.org = syn2ann_thp(nss.org);
sps.org = syn2ann_thp(sps.org);
nss.org = syn2ann_spp(nss.org);
sps.org = syn2ann_spp(sps.org);

%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Hybridization\n');
[nss.hyb,sps.hyb,hbs_org] = lfhf_mashup(nss.org,sps.org);
% hbs = syn2ann_blc(hbs_org);
hbs = hbs_org;

close all;
figure
loglog(BB_FREQ,abs(BB_AMP)); hold all;
loglog(hbs.mon.vfr{1},abs(hbs.syn{1}.fsa.e),'r--'); 
loglog(nss.hyb.mon.vfr{1},abs(nss.hyb.syn{1}.fsa.e),'c--');format_figures;

figure
plot(BB_TIME,BB_ACC); hold all;
plot(hbs.mon.vtm{1},hbs.syn{1}.tha.e,'r--'); format_figures;

figure
plot(BB_TIME,BB_VEL); hold all;
plot(hbs.mon.vtm{1},hbs.syn{1}.thv.e,'r--'); format_figures;

figure
plot(BB_TIME,BB_DIS); hold all;
plot(hbs.mon.vtm{1},hbs.syn{1}.thd.e,'r--'); format_figures;

figure
plot(nss.org.mon.vtm{1},nss.org.syn{1}.tha.e);hold all;
plot(nss.hyb.mon.vtm{1},nss.hyb.syn{1}.tha.e);hold all;
format_figures;

figure
plot(nss.org.mon.vtm{1},nss.org.syn{1}.thv.e);hold all;
plot(nss.hyb.mon.vtm{1},nss.hyb.syn{1}.thv.e);hold all;
format_figures;

figure
plot(nss.org.mon.vtm{1},nss.org.syn{1}.thd.e);hold all;
plot(nss.hyb.mon.vtm{1},nss.hyb.syn{1}.thd.e);hold all;
format_figures;

figure
plot(sps.org.mon.vtm{1},sps.org.syn{1}.tha.e);hold all;
plot(sps.hyb.mon.vtm{1},sps.hyb.syn{1}.tha.e);hold all;
format_figures;

figure
plot(sps.org.mon.vtm{1},sps.org.syn{1}.thv.e);hold all;
plot(sps.hyb.mon.vtm{1},sps.hyb.syn{1}.thv.e);hold all;
format_figures;

figure
plot(sps.org.mon.vtm{1},sps.org.syn{1}.thd.e);hold all;
plot(sps.hyb.mon.vtm{1},sps.hyb.syn{1}.thd.e);hold all;
format_figures;

keyboard

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
nss.hyb = syn2ann_thp(nss.hyb);
sps.hyb = syn2ann_thp(sps.hyb);
hbs_org = syn2ann_thp(hbs_org);
hbs     = syn2ann_thp(hbs);

%% *SPECTRA*
fprintf('--> Spectra\n');
nss.hyb = syn2ann_spp(nss.hyb);
sps.hyb = syn2ann_spp(sps.hyb);
hbs_org = syn2ann_spp(hbs_org);
hbs = syn2ann_spp(hbs);

