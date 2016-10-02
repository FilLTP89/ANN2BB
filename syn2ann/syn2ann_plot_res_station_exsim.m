global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
fprintf('%s-with HYBRID\n',bhr.nm{mm_});
dlg = cell(3,1);
dlg{1} = 'EW';
dlg{2} = 'NS';
dlg{3} = 'UD';
    
for j_ = 1:numel(cpp)
    %% *POST-PROCESS - NUMERICAL SIMULATIONS VS RECORDS*
    
    fn = fullfile(sp,sprintf('%s_sim_rec_%s',st,ann.mtd.scl));
    % all
    set(0,'defaultaxescolororder',clr0);
    syn2ann_plot_compare(...
        ones(5,1),...
        {rec.org,nss.org},...
        mm_,dlg(j_),...
        {'REC';'PBS'},...
        (cpp{j_}),...
        fn,...
        [0,vtm_shift(2)]);
    
    %% *POST-PROCESS - BB HYBRIDIZATION*
    
    fn = fullfile(sp,sprintf('%s_sim_exsim_hyb_%s',st,ann.mtd.scl));
    % time-histories
    set(0,'defaultaxescolororder',clr10);
    syn2ann_plot_compare(...
        [0 0 1 1 1],...
        {nss.hyb.exs;exs.hyb;hbs.exs},...
        mm_,dlg(j_),...
        {'PBS';'EXSIM';'HYB'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(2)*ones(3,1));
    % fourier spectra
    set(0,'defaultaxescolororder',clr11);
    syn2ann_plot_compare(...
        [0 1 0 0 0],...
        {nss.exs;exs.org;nss.hyb.exs;exs.hyb;hbs.exs},...
        mm_,dlg(j_),...
        {'';'';'PBS';'EXSIM';'HYB'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(2)*ones(5,1),...
        {'none','none','none','none','none'},...
        {'--','--','-','-','-'});
    % psa
    set(0,'defaultaxescolororder',clr12);
    syn2ann_plot_compare(...
        [1 0 0 0 0],...
        {rec.org;nss.hyb.exs;exs.hyb;hbs.exs},...
        mm_,dlg(j_),...
        {'REC';'PBS';'EXSIM';'HYB'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(2)*ones(5,1),...
        {'none','none','none','none','none'},...
        {'-','--','--','-'});
    
    %% *POST-PROCESS - HYBRIDS vs ANN vs SPECTRAL MATCHED*
    
    fn = fullfile(sp,sprintf('%s_hyb_ann_exsim_%s',st,ann.mtd.scl));
    set(0,'defaultaxescolororder',clr2);
    syn2ann_plot_compare(...
        [1 0 0 0 0],...
        {rec.org;hbs.exs;trs.exs.(cpp{j_});spm.exs.(cpp{j_})},...
        mm_,dlg(j_),...
        {'REC';'HYB';'ANN';'SPM'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(1)*ones(4,1),...
        {'none';'none';'o';'none'},...
        {'-';'-';'none';'-'});
    
    %% *POST-PROCESS - SPECTRAL MATCHING vs RECORDS*
    
    fn = fullfile(sp,sprintf('%s_spm_rec_exsim_%s',st,ann.mtd.scl));
    set(0,'defaultaxescolororder',clr3);
    syn2ann_plot_compare(...
        [0 1 1 1 1],...
        {rec.org;spm.exs.(cpp{j_})},...
        mm_,dlg(j_),...
        {'REC';'SPM'},...
        (cpp{j_}),...
        fn,...
        [0,vtm_shift(1)]);
end