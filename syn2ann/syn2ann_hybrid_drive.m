%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');
% % ============ DEBUG CODE---> CHECK WITH BBSYNT ============================
% [BB_DTIME,BB_ACC,BB_FREQ,BB_AMP] = bbsynt_original(pbs.org.mon.dtm(1),pbs.org.syn{1}.tha.e,...
%     sps.org.mon.dtm(1),sps.org.syn{1}.tha.e,mon.fa(1),mon.fb(1));
% BB_TIME = BB_DTIME*(0:numel(BB_ACC)-1);
% BB_VEL = cumtrapz(BB_ACC)*BB_DTIME;
% BB_DIS = cumtrapz(BB_VEL)*BB_DTIME;
% BB_FREQ(end+1:numel(BB_AMP)) = BB_FREQ(end)+mean(diff(BB_FREQ))*(1:numel(BB_AMP)/2)';

%
% _SABETTA & PUGLIESE 1996_
%
%% *RESAMPLING*
fprintf('--> Resampling\n');
[pbs.sps,sps.org] = lfhf_rsmpl(pbs.org,sps.org);

%% *PAD LF/HF*
fprintf('--> Padding/Tapering\n');
[pbs.sps,sps.org] = lfhf_pad(pbs.sps,sps.org);

%% *ALIGN LF/HF*
fprintf('--> Align records\n');
[pbs.sps,sps.org] = lfhf_shift(pbs.sps,sps.org);
pbs.sps = syn2ann_thp(pbs.sps);
sps.org = syn2ann_thp(sps.org);
pbs.sps = syn2ann_spp(pbs.sps);
sps.org = syn2ann_spp(sps.org);

%% *SPECTRAL MASHUP LF/HF*
fprintf('--> Hybridization\n');
[pbs.hyb.sps,sps.hyb,hbs.sps] = lfhf_mashup(pbs.sps,sps.org);

%% *PGA-PGV-PGD & ARIAS INTENSITY*
fprintf('--> Peak Values and Arias\n');
pbs.hyb.sps = syn2ann_thp(pbs.hyb.sps);
sps.hyb = syn2ann_thp(sps.hyb);
hbs.sps = syn2ann_thp(hbs.sps);

%% *SPECTRA*
fprintf('--> Spectra\n');
pbs.hyb.sps = syn2ann_spp(pbs.hyb.sps);
sps.hyb = syn2ann_spp(sps.hyb);
hbs.sps = syn2ann_spp(hbs.sps);
% switch lower(hybrid_type)
%     case 'sp96'
%         %
%         % _SABETTA & PUGLIESE 1996_
%         %
%         %% *RESAMPLING*
%         fprintf('--> Resampling\n');
%         [pbs.sps,sps.org] = lfhf_rsmpl(pbs.org,sps.org);
%         
%         %% *PAD LF/HF*
%         fprintf('--> Padding/Tapering\n');
%         [pbs.sps,sps.org] = lfhf_pad(pbs.sps,sps.org);
%         
%         %% *ALIGN LF/HF*
%         fprintf('--> Align records\n');
%         [pbs.sps,sps.org] = lfhf_shift(pbs.sps,sps.org);
%         pbs.sps = syn2ann_thp(pbs.sps);
%         sps.org = syn2ann_thp(sps.org);
%         pbs.sps = syn2ann_spp(pbs.sps);
%         sps.org = syn2ann_spp(sps.org);
%         
%         %% *SPECTRAL MASHUP LF/HF*
%         fprintf('--> Hybridization\n');
%         [pbs.hyb.sps,sps.hyb,hbs.sps] = lfhf_mashup(pbs.sps,sps.org);
%         
%         %% *PGA-PGV-PGD & ARIAS INTENSITY*
%         fprintf('--> Peak Values and Arias\n');
%         pbs.hyb.sps = syn2ann_thp(pbs.hyb.sps);
%         sps.hyb = syn2ann_thp(sps.hyb);
%         hbs.sps = syn2ann_thp(hbs.sps);
%         
%         %% *SPECTRA*
%         fprintf('--> Spectra\n');
%         pbs.hyb.sps = syn2ann_spp(pbs.hyb.sps);
%         sps.hyb = syn2ann_spp(sps.hyb);
%         hbs.sps = syn2ann_spp(hbs.sps);
%         
%     case 'exsim'
%         %
%         % _EXSIM_
%         %
%         [pbs.exs,exs.org] = lfhf_rsmpl(pbs.org,exs.org);
%         
%         %% *PAD LF/HF*
%         fprintf('--> Padding/Tapering\n');
%         [pbs.exs,exs.org] = lfhf_pad(pbs.exs,exs.org);
%         
%         %% *ALIGN LF/HF*
%         fprintf('--> Align records\n');
%         [pbs.exs,exs.org] = lfhf_shift(pbs.exs,exs.org);
%         pbs.exs = syn2ann_thp(pbs.exs);
%         exs.org = syn2ann_thp(exs.org);
%         pbs.exs = syn2ann_spp(pbs.exs);
%         exs.org = syn2ann_spp(exs.org);
%         
%         %% *SPECTRAL MASHUP LF/HF*
%         fprintf('--> Hybridization\n');
%         [pbs.hyb.exs,exs.hyb,hbs.exs] = lfhf_mashup(pbs.exs,exs.org);
%         
%         %% *PGA-PGV-PGD & ARIAS INTENSITY*
%         fprintf('--> Peak Values and Arias\n');
%         pbs.hyb.exs = syn2ann_thp(pbs.hyb.exs);
%         exs.hyb = syn2ann_thp(exs.hyb);
%         hbs.exs = syn2ann_thp(hbs.exs);
%         
%         %% *SPECTRA*
%         fprintf('--> Spectra\n');
%         pbs.hyb.exs = syn2ann_spp(pbs.hyb.exs);
%         exs.hyb = syn2ann_spp(exs.hyb);
%         hbs.exs = syn2ann_spp(hbs.exs);
%         
%     case 'both'
%         %% *RESAMPLING*
%         fprintf('--> Resampling\n');
%         [pbs.sps,sps.org] = lfhf_rsmpl(pbs.org,sps.org);
%         [pbs.exs,exs.org] = lfhf_rsmpl(pbs.org,exs.org);
%         %% *PAD LF/HF*
%         fprintf('--> Padding/Tapering\n');
%         [pbs.sps,sps.org] = lfhf_pad(pbs.sps,sps.org);
%         [pbs.exs,exs.org] = lfhf_pad(pbs.exs,exs.org);
%         %% *ALIGN LF/HF*
%         fprintf('--> Align records\n');
%         [pbs.sps,sps.org] = lfhf_shift(pbs.sps,sps.org);
%         pbs.sps = syn2ann_thp(pbs.sps);
%         sps.org = syn2ann_thp(sps.org);
%         pbs.sps = syn2ann_spp(pbs.sps);
%         sps.org = syn2ann_spp(sps.org);
%         %
%         [pbs.exs,exs.org] = lfhf_shift(pbs.exs,exs.org);
%         pbs.exs = syn2ann_thp(pbs.exs);
%         exs.org = syn2ann_thp(exs.org);
%         pbs.exs = syn2ann_spp(pbs.exs);
%         exs.org = syn2ann_spp(exs.org);
%         
%         %% *SPECTRAL MASHUP LF/HF*
%         fprintf('--> Hybridization\n');
%         [pbs.hyb.sps,sps.hyb,hbs.sps] = lfhf_mashup(pbs.sps,sps.org);
%         [pbs.hyb.exs,exs.hyb,hbs.exs] = lfhf_mashup(pbs.exs,exs.org);
%         
%         %% *PGA-PGV-PGD & ARIAS INTENSITY*
%         fprintf('--> Peak Values and Arias\n');
%         pbs.hyb.sps = syn2ann_thp(pbs.hyb.sps);
%         sps.hyb = syn2ann_thp(sps.hyb);
%         hbs.sps = syn2ann_thp(hbs.sps);
%         
%         %% *SPECTRA*
%         fprintf('--> Spectra\n');
%         pbs.hyb.sps = syn2ann_spp(pbs.hyb.sps);
%         sps.hyb = syn2ann_spp(sps.hyb);
%         hbs.sps = syn2ann_spp(hbs.sps);        
%         
%         %% *PGA-PGV-PGD & ARIAS INTENSITY*
%         fprintf('--> Peak Values and Arias\n');
%         pbs.hyb.exs = syn2ann_thp(pbs.hyb.exs);
%         exs.hyb = syn2ann_thp(exs.hyb);
%         hbs.exs = syn2ann_thp(hbs.exs);
%         
%         %% *SPECTRA*
%         fprintf('--> Spectra\n');
%         pbs.hyb.exs = syn2ann_spp(pbs.hyb.exs);
%         exs.hyb = syn2ann_spp(exs.hyb);
%         hbs.exs = syn2ann_spp(hbs.exs);
% end



% close all;
% figure
% loglog(BB_FREQ,abs(BB_AMP)); hold all;
% loglog(hbs.mon.vfr{1},abs(hbs.syn{1}.fsa.e),'r--');
% loglog(pbs.hyb.mon.vfr{1},abs(pbs.hyb.syn{1}.fsa.e),'c--');format_figures;
%
% figure
% plot(BB_TIME,BB_ACC); hold all;
% plot(hbs.mon.vtm{1},hbs.syn{1}.tha.e,'r--'); format_figures;
%
% figure
% plot(BB_TIME,BB_VEL); hold all;
% plot(hbs.mon.vtm{1},hbs.syn{1}.thv.e,'r--'); format_figures;
%
% figure
% plot(BB_TIME,BB_DIS); hold all;
% plot(hbs.mon.vtm{1},hbs.syn{1}.thd.e,'r--'); format_figures;
%
% figure
% plot(pbs.org.mon.vtm{1},pbs.org.syn{1}.tha.e);hold all;
% plot(pbs.hyb.mon.vtm{1},pbs.hyb.syn{1}.tha.e);hold all;
% format_figures;
%
% figure
% plot(pbs.org.mon.vtm{1},pbs.org.syn{1}.thv.e);hold all;
% plot(pbs.hyb.mon.vtm{1},pbs.hyb.syn{1}.thv.e);hold all;
% format_figures;
%
% figure
% plot(pbs.org.mon.vtm{1},pbs.org.syn{1}.thd.e);hold all;
% plot(pbs.hyb.mon.vtm{1},pbs.hyb.syn{1}.thd.e);hold all;
% format_figures;
%
% figure
% plot(sps.org.mon.vtm{1},sps.org.syn{1}.tha.e);hold all;
% plot(sps.hyb.mon.vtm{1},sps.hyb.syn{1}.tha.e);hold all;
% format_figures;
%
% figure
% plot(sps.org.mon.vtm{1},sps.org.syn{1}.thv.e);hold all;
% plot(sps.hyb.mon.vtm{1},sps.hyb.syn{1}.thv.e);hold all;
% format_figures;
%
% figure
% plot(sps.org.mon.vtm{1},sps.org.syn{1}.thd.e);hold all;
% plot(sps.hyb.mon.vtm{1},sps.hyb.syn{1}.thd.e);hold all;
% format_figures;
%
% keyboard


