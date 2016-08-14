%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');
% ============ DEBUG CODE---> CHECK WITH BBSYNT ============================
[BB_DTIME,BB_ACC,BB_FREQ,BB_AMP] = bbsynt_original(nss.org.mon.dtm(1),nss.org.syn{1}.tha.e,...
    sps.org.mon.dtm(1),sps.org.syn{1}.tha.e,mon.fa(1),mon.fb(1));
BB_TIME = BB_DTIME*(0:numel(BB_ACC)-1);
disp('==== debug ======')
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
[nss.hyb,sps.hyb,hbs] = lfhf_mashup(nss.org,sps.org);


figure
BB_FREQ(end+1:numel(BB_AMP))=BB_FREQ(end)+mean(diff(BB_FREQ))*(1:numel(BB_AMP)/2)';
loglog(BB_FREQ,abs(BB_AMP)); hold all;
loglog(hbs.mon.vfr{1},abs(hbs.syn{1}.fsa.e),'r--'); 
loglog(nss.hyb.mon.vfr{1},abs(nss.hyb.syn{1}.fsa.e),'c--');format_figures;
xlim([50,100]);
% figure
% plot(BB_TIME,BB_ACC); hold all;
% plot(hbs.mon.vtm{1},hbs.syn{1}.tha.e,'r--'); format_figures;
% 
keyboard

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
nss.hyb = syn2ann_thp(nss.hyb);
sps.hyb = syn2ann_thp(sps.hyb);
hbs     = syn2ann_thp(hbs);

%% *SPECTRA*
fprintf('--> Spectra\n');
nss.hyb = syn2ann_spp(nss.hyb);
sps.hyb = syn2ann_spp(sps.hyb);
hbs = syn2ann_spp(hbs);
