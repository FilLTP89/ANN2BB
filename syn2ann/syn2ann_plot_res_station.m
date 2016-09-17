global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
fprintf('%s-with HYBRID\n',bhr.nm{mm_});
dlg = cell(3,1);
dlg{1} = 'EW';
dlg{2} = 'NS';
dlg{3} = 'UD';
    
for j_ = 1:numel(cpp)
    %% *POST-PROCESS - ORIGINAL VS FILTERED RECORDS*
    fn = fullfile(sp,sprintf('%s_rec_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        ones(5,1),...
        {rec.org},...
        mm_,dlg(j_),...
        {'RECORDED'},...
        (cpp{j_}),...
        fn,...
        0);
    %% *POST-PROCESS - NUMERICAL SIMULATIONS VS RECORDS*
    fn = fullfile(sp,sprintf('%s_sim_rec_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        ones(5,1),...
        {rec.org,nss.org},...
        mm_,dlg(j_),...
        {'RECORDED';'SIMULATED'},...
        (cpp{j_}),...
        fn,...
        [0,vtm_shift(2)]);
    
    %% *POST-PROCESS - BB HYBRIDIZATION*
    fn = fullfile(sp,sprintf('%s_sim_sp96_hyb_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        [1 0 1 1 1],...
        {nss.hyb;sps.hyb;hbs_org},...
        mm_,dlg(j_),...
        {'SIMULATED';'SP96';'HYBRID'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(2)*ones(3,1));
    set(0,'defaultaxescolororder',clr1);
    syn2ann_plot_compare(...
        [0 1 0 0 0],...
        {nss.org;sps.org;nss.hyb;sps.hyb;hbs_org},...
        mm_,dlg(j_),...
        {'';'';'SIMULATED';'SP96';'HYBRID'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(2)*ones(5,1),...
        {'none','none','none','none','none'},...
        {'--','--','-','-','-'});
    set(0,'defaultaxescolororder',clr);
    %% *POST-PROCESS - HYBRIDS vs ANN*
    fn = fullfile(sp,sprintf('%s_hyb_ann_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        [1 0 0 0 0],...
        {hbs;trs.(cpp{j_})},...
        mm_,dlg(j_),...
        {'HYBRID';'ANN'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(1)*ones(2,1),...
        {'none';'o'},...
        {'-';'none'});
    
    %% *POST-PROCESS - SPECTRAL MATCHING vs ANN*
    fn = fullfile(sp,sprintf('%s_spm_ann_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        [1 0 0 0 0],...
        {spm.(cpp{j_});...
        trs.(cpp{j_})},...
        mm_,dlg(j_),...
        {'MATCHED';'ANN'},...
        (cpp{j_}),...
        fn,...
        vtm_shift(1)*ones(2,1),...
        {'none';'o'},...
        {'-';'none'});
    %
    %% *POST-PROCESS - SPECTRAL MATCHING vs RECORDS*
    fn = fullfile(sp,sprintf('%s_spm_rec_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        [1 0 0 1 1],...
        {rec.fil;spm.(cpp{j_})},...
        mm_,dlg(j_),...
        {'RECORDED';'MATCHED'},...
        (cpp{j_}),...
        fn,...
        [0,vtm_shift(1)]);
    
    %% *POST-PROCESS - SPECTRAL MATCHING vs RECORDS*
    fn = fullfile(sp,sprintf('%s_spm_rec_%s',st,ann.mtd.scl));
    syn2ann_plot_compare(...
        [0 1 1 0 0],...
        {rec.fil;spm.(cpp{j_})},...
        mm_,dlg(j_),...
        {'RECORDED';'MATCHED'},...
        (cpp{j_}),...
        fn,...
        [0,vtm_shift(1)]);
end