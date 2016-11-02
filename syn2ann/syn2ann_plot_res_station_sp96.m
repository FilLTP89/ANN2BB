global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
fprintf('%s-with HYBRID\n',bhr.nm{mm_});
dlg = cell(3,1);
dlg{1} = 'EW';
dlg{2} = 'NS';
dlg{3} = 'UD';
    
for j_ = 1:numel(cpp)
    %% *POST-PROCESS - RECORDS vs PBS*
    fnm = fullfile(sp,sprintf('%s_rec_pbs_%s',st,ann.mtd.scl{j_}));
    % all
    set(0,'defaultaxescolororder',clr.rec_pbs);
    syn2ann_plot_compare(...
        ones(5,1),...
        {rec.org,nss.org},...
        mm_,dlg(j_),...
        {'REC';'PBS'},...
        (cpp{j_}),bhr.nm{mm_},...
        fnm,...
        [0,vtm_shift(2)]);
    
    %% *POST-PROCESS - BB HYBRIDIZATION (EMP)*
    fnm = fullfile(sp,sprintf('%s_pbs_emp_hyb_%s',st,ann.mtd.scl{j_}));
    % fourier spectra
    set(0,'defaultaxescolororder',clr.pbs_emp_hyb_fss);
    syn2ann_plot_compare(...
        [0 1 0 0 0],...
        {nss.sps;sps.org;nss.hyb.sps;sps.hyb;hbs.sps},...
        mm_,dlg(j_),...
        {'';'';'PBS';'SP96';'HYB'},...
        (cpp{j_}),bhr.nm{mm_},...
        fnm,...
        vtm_shift(2)*ones(5,1),...
        {'none','none','none','none','none'},...
        {':',':','-','-','-'});
    % psa
    set(0,'defaultaxescolororder',clr.rec_pbs_emp_hyb);
    syn2ann_plot_compare(...
        [1 0 0 0 0],...
        {rec.org;nss.hyb.sps;sps.hyb;hbs.sps},...
        mm_,dlg(j_),...
        {'REC';'PBS';'SP96';'HYB'},...
        (cpp{j_}),bhr.nm{mm_},...
        fnm,...
        vtm_shift(2)*ones(5,1),...
        {'none','none','none','none','none'},...
        {'-',':','-.','-'});
    
    %% *POST-PROCESS - RECORDS vs SPECTRAL MATCHED (EMP)*
    fnm = fullfile(sp,sprintf('%s_rec_spm_EMP_%s',st,ann.mtd.scl{j_}));
    set(0,'defaultaxescolororder',clr.rec_spm);
    syn2ann_plot_compare(...
        [0 0 1 1 1],...
        {rec.org;spm.sps.(cpp{j_})},...
        mm_,dlg(j_),...
        {'REC';'SPM'},...
        (cpp{j_}),bhr.nm{mm_},...
        fnm,...
        [0,vtm_shift(1)]);
    %% *POST-PROCESS - RECORDS vs PBS vs SPECTRAL MATCHED (EMP)*
    fnm = fullfile(sp,sprintf('%s_rec_pbs_spm_EMP_%s',st,ann.mtd.scl{j_}));
    set(0,'defaultaxescolororder',clr.rec_pbs_spm);
    syn2ann_plot_compare(...
        [1 1 0 0 0],...
        {rec.org;nss.sps;spm.sps.(cpp{j_})},...
        mm_,dlg(j_),...
        {'REC';'PBS';'SPM'},...
        (cpp{j_}),bhr.nm{mm_},...
        fnm,[0,vtm_shift(1),vtm_shift(1)],...
        {'none';'none';'none'},...
        {'-',':','-'});
end