%% *LF-HF HYBRIDIZATION*
fprintf('---------------------\n4. HYBRIDIZATION\n---------------------\n');
% % ============ DEBUG CODE---> CHECK WITH BBSYNT ============================
% [BB_DTIME,BB_ACC,BB_FREQ,BB_AMP] = bbsynt_original(pbs.org.mon.dtm(1),pbs.org.syn{1}.tha.e,...
%     sps.org.mon.dtm(1),sps.org.syn{1}.tha.e,mon.fa(1),mon.fb(1));
% BB_TIME = BB_DTIME*(0:numel(BB_ACC)-1);
% BB_VEL = cumtrapz(BB_ACC)*BB_DTIME;
% BB_DIS = cumtrapz(BB_VEL)*BB_DTIME;
% BB_FREQ(end+1:numel(BB_AMP)) = BB_FREQ(end)+mean(diff(BB_FREQ))*(1:numel(BB_AMP)/2)';

switch lower(hybrid_type)
    case 'sp96'
        %
        % _SABETTA & PUGLIESE 1996_
        %
        %% *RESAMPLING*
        fprintf('--> Resampling\n');
        [nss.sps,sps.org] = lfhf_rsmpl(nss.org,sps.org);
        
        %% *PAD LF/HF*
        fprintf('--> Padding/Tapering\n');
        [nss.sps,sps.org] = lfhf_pad(nss.sps,sps.org);
        
        %% *ALIGN LF/HF*
        fprintf('--> Align records\n');
        [nss.sps,sps.org] = lfhf_shift(nss.sps,sps.org);
        nss.sps = syn2ann_thp(nss.sps);
        sps.org = syn2ann_thp(sps.org);
        nss.sps = syn2ann_spp(nss.sps);
        sps.org = syn2ann_spp(sps.org);
        
        %% *SPECTRAL MASHUP LF/HF*
        fprintf('--> Hybridization\n');
        [nss.hyb.sps,sps.hyb,hbs.sps] = lfhf_mashup(nss.sps,sps.org);
        
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        nss.hyb.sps = syn2ann_thp(nss.hyb.sps);
        sps.hyb = syn2ann_thp(sps.hyb);
        hbs.sps = syn2ann_thp(hbs.sps);
        
        %% *SPECTRA*
        fprintf('--> Spectra\n');
        nss.hyb.sps = syn2ann_spp(nss.hyb.sps);
        sps.hyb = syn2ann_spp(sps.hyb);
        hbs.sps = syn2ann_spp(hbs.sps);
        
    case 'exsim'
        %
        % _EXSIM_
        %
        [nss.exs,exs.org] = lfhf_rsmpl(nss.org,exs.org);
        
        %% *PAD LF/HF*
        fprintf('--> Padding/Tapering\n');
        [nss.exs,exs.org] = lfhf_pad(nss.exs,exs.org);
        
        %% *ALIGN LF/HF*
        fprintf('--> Align records\n');
        [nss.exs,exs.org] = lfhf_shift(nss.exs,exs.org);
        nss.exs = syn2ann_thp(nss.exs);
        exs.org = syn2ann_thp(exs.org);
        nss.exs = syn2ann_spp(nss.exs);
        exs.org = syn2ann_spp(exs.org);
        
        %% *SPECTRAL MASHUP LF/HF*
        fprintf('--> Hybridization\n');
        [nss.hyb.exs,exs.hyb,hbs.exs] = lfhf_mashup(nss.exs,exs.org);
        
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        nss.hyb.exs = syn2ann_thp(nss.hyb.exs);
        exs.hyb = syn2ann_thp(exs.hyb);
        hbs.exs = syn2ann_thp(hbs.exs);
        
        %% *SPECTRA*
        fprintf('--> Spectra\n');
        nss.hyb.exs = syn2ann_spp(nss.hyb.exs);
        exs.hyb = syn2ann_spp(exs.hyb);
        hbs.exs = syn2ann_spp(hbs.exs);
        
    case 'both'
        %% *RESAMPLING*
        fprintf('--> Resampling\n');
        [nss.sps,sps.org] = lfhf_rsmpl(nss.org,sps.org);
        [nss.exs,exs.org] = lfhf_rsmpl(nss.org,exs.org);
        %% *PAD LF/HF*
        fprintf('--> Padding/Tapering\n');
        [nss.sps,sps.org] = lfhf_pad(nss.sps,sps.org);
        [nss.exs,exs.org] = lfhf_pad(nss.exs,exs.org);
        %% *ALIGN LF/HF*
        fprintf('--> Align records\n');
        [nss.sps,sps.org] = lfhf_shift(nss.sps,sps.org);
        nss.sps = syn2ann_thp(nss.sps);
        sps.org = syn2ann_thp(sps.org);
        nss.sps = syn2ann_spp(nss.sps);
        sps.org = syn2ann_spp(sps.org);
        %
        [nss.exs,exs.org] = lfhf_shift(nss.exs,exs.org);
        nss.exs = syn2ann_thp(nss.exs);
        exs.org = syn2ann_thp(exs.org);
        nss.exs = syn2ann_spp(nss.exs);
        exs.org = syn2ann_spp(exs.org);
        
        %% *SPECTRAL MASHUP LF/HF*
        fprintf('--> Hybridization\n');
        [nss.hyb.sps,sps.hyb,hbs.sps] = lfhf_mashup(nss.sps,sps.org);
        [nss.hyb.exs,exs.hyb,hbs.exs] = lfhf_mashup(nss.exs,exs.org);
        
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        nss.hyb.sps = syn2ann_thp(nss.hyb.sps);
        sps.hyb = syn2ann_thp(sps.hyb);
        hbs.sps = syn2ann_thp(hbs.sps);
        
        %% *SPECTRA*
        fprintf('--> Spectra\n');
        nss.hyb.sps = syn2ann_spp(nss.hyb.sps);
        sps.hyb = syn2ann_spp(sps.hyb);
        hbs.sps = syn2ann_spp(hbs.sps);        
        
        %% *PGA-PGV-PGD & ARIAS INTENSITY*
        fprintf('--> Peak Values and Arias\n');
        nss.hyb.exs = syn2ann_thp(nss.hyb.exs);
        exs.hyb = syn2ann_thp(exs.hyb);
        hbs.exs = syn2ann_thp(hbs.exs);
        
        %% *SPECTRA*
        fprintf('--> Spectra\n');
        nss.hyb.exs = syn2ann_spp(nss.hyb.exs);
        exs.hyb = syn2ann_spp(exs.hyb);
        hbs.exs = syn2ann_spp(hbs.exs);
end



% close all;
% figure
% loglog(BB_FREQ,abs(BB_AMP)); hold all;
% loglog(hbs.mon.vfr{1},abs(hbs.syn{1}.fsa.e),'r--');
% loglog(nss.hyb.mon.vfr{1},abs(nss.hyb.syn{1}.fsa.e),'c--');format_figures;
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
% plot(nss.org.mon.vtm{1},nss.org.syn{1}.tha.e);hold all;
% plot(nss.hyb.mon.vtm{1},nss.hyb.syn{1}.tha.e);hold all;
% format_figures;
%
% figure
% plot(nss.org.mon.vtm{1},nss.org.syn{1}.thv.e);hold all;
% plot(nss.hyb.mon.vtm{1},nss.hyb.syn{1}.thv.e);hold all;
% format_figures;
%
% figure
% plot(nss.org.mon.vtm{1},nss.org.syn{1}.thd.e);hold all;
% plot(nss.hyb.mon.vtm{1},nss.hyb.syn{1}.thd.e);hold all;
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


