%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');

%% *RESAMPLING*
fprintf('--> Resampling\n');
[nss.org,sps.org] = lfhf_rsmpl(nss.org,sps.org);

%% *PAD LF/HF*
fprintf('--> Padding/Tapering\n');
[nss.org,sps.org] = lfhf_pad(nss.org,sps.org);

%% *ALIGN LF/HF*
fprintf('--> Align records\n');
[nss.org,sps.org] = lfhf_shift(nss.org,sps.org);
% nss.org = syn2ann_thp(nss.org);
% sps.org = syn2ann_thp(sps.org);
nss.org = syn2ann_spp(nss.org,2);
sps.org = syn2ann_spp(sps.org,2);

%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Hybridization\n');
[nss.hyb,sps.hyb,hbs] = lfhf_mashup(nss.org,sps.org);

% %% *PGA-PGV-PGD & ARIAS INTENSITY*
% fprintf('--> Peak Values and Arias\n');
% nss.hyb = syn2ann_thp(nss.hyb);
% sps.hyb = syn2ann_thp(sps.hyb);
hbs     = syn2ann_thp(hbs);

% %% *SPECTRA*
% fprintf('--> Spectra\n');
% nss.hyb = syn2ann_spp(nss.hyb);
% sps.hyb = syn2ann_spp(sps.hyb);
hbs = syn2ann_spp(hbs);

%============ DEBUG CODE---> CHECK WITH BBSYNT ============================
% [BB_DTIME,BB_ACC,BB_FREQ,BB_AMP] = bbsynt_original(nss.org.mon.dtm(1),nss.org.syn{1}.tha.e,...
%     sps.org.mon.dtm(1),sps.org.syn{1}.tha.e,mon.fa,mon.fb);
% BB_TIME = BB_DTIME*(0:numel(BB_ACC)-1);
% disp('==== debug ======')
% figure
% loglog(BB_FREQ,abs(BB_AMP(1:numel(BB_FREQ)))); hold all;
% loglog(hbs.mon.vfr{1},abs(hbs.syn{1}.fsa.e),'r--'); format_figures;
% 
% figure
% plot(BB_TIME,BB_ACC); hold all;
% plot(hbs.mon.vtm{1},hbs.syn{1}.tha.e,'r--'); format_figures;
% keyboard