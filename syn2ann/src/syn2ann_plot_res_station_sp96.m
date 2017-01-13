global pfg xlm xlb xtk ylm ylb ytk grd scl tit utd 

for j_ = 1:numel(cpp)
    %% *POST-PROCESS - RECORDS vs PBS*
    set(0,'defaultaxescolororder',rgb('Black'));
    fnm = fullfile(sp,sprintf('%s_rec_%s',st,ann.(cpp{j_}).nl));
    syn2ann_plot_compare([0;0;1;1;1],{rec.org},mm_,dlg(j_),...
        {'REC'},(cpp{j_}),ttt,fnm,0,{'none'},...
        {'-'},2);
    fnm = fullfile(sp,sprintf('%s_pbs_%s',st,ann.(cpp{j_}).nl));
    syn2ann_plot_compare([0;0;1;1;1],{pbs.org},mm_,dlg(j_),...
        {'PBS'},(cpp{j_}),{' '},fnm,vtm_shift,{'none'},...
        {'-'},2);
    fnm = fullfile(sp,sprintf('%s_hsp_%s',st,ann.(cpp{j_}).nl));
    syn2ann_plot_compare([0;0;1;1;1],{hbs.bst},mm_,dlg(j_),...
        {'HYB-SP96'},(cpp{j_}),{' '},fnm,vtm_shift,{'none'},...
        {'-'},2);
    fnm = fullfile(sp,sprintf('%s_ssp_%s',st,ann.(cpp{j_}).nl));
    syn2ann_plot_compare([0;0;1;1;1],{spm.sps.(cpp{j_})},mm_,dlg(j_),...
        {'SPM'},(cpp{j_}),{' '},fnm,vtm_shift,{'none'},...
        {'-'},2);
    set(0,'defaultaxescolororder',clr131);
    fnm = fullfile(sp,sprintf('%s_rec_pbs_hsp_ssp_%s',st,ann.(cpp{j_}).nl));
    syn2ann_plot_compare([1;0;0;0;0],{rec.org;pbs.org;hbs.bst;spm.sps.(cpp{j_})},...
        mm_,dlg(j_),{'REC';'PBS';'HYB';'SPM'},(cpp{j_}),ttt,fnm,zeros(1,4),...
        {'none';'none';'none';'none'},...
        {'-';':';'--';'-'},...
        [4,1.5,1,2.5]);
    fnm = fullfile(sp,sprintf('%s_rec_pbs_hsp_ssp_%s',st,ann.(cpp{j_}).nl));
    syn2ann_plot_compare([0;1;0;0;0],{rec.org;pbs.org;hbs.bst;spm.sps.(cpp{j_})},...
        mm_,dlg(j_),{'REC';'PBS';'HYB';'SPM'},(cpp{j_}),ttt,fnm,zeros(1,4),...
        {'none';'none';'none';'none'},...
        {'-';':';'--';'-'},...
        [1,1,1,1]);
%     %% *POST-PROCESS - RECORDS vs HYB vs SPM*
%     set(0,'defaultaxescolororder',rgb('Black'));
%     fnm = fullfile(sp,sprintf('%s_rec_hyb_spm_%s',st,ann.mtd.scl{j_}));
%     syn2ann_plot_compare([0;0;1;0;0],{rec.org;hbs.bst;spm.sps.(cpp{j_})},mm_,dlg(j_),...
%         {'REC';'PBS';'SPM'},(cpp{j_}),ttt,fnm,[0;vtm_shift;vtm_shift],{'none';'none',;'none'},...
%         {'-';'-';'-'},[2,2,2]);
%     % PBS-TIME HISTORIES
%     fnm = fullfile(sp,sprintf('%s_pbs_%s',st,ann.mtd.scl{j_}));
%     syn2ann_plot_compare([0;0;1;1;1],{pbs.org},mm_,dlg(j_),{'PBS'},(cpp{j_}),...
%         bhr.nm{mm_},fnm,vtm_shift,{'none'},{'-'},2);
    
%     %% *POST-PROCESS - BB HYBRIDIZATION (EMP)*
%     fnm = fullfile(sp,sprintf('%s_pbs_emp_hyb_%s',st,ann.mtd.scl{j_}));
%     % fourier spectra
%     set(0,'defaultaxescolororder',clr.pbs_emp_hyb_fss);
%     syn2ann_plot_compare(...
%         [0 1 0 0 0],...
%         {pbs.bst;sps.org;pbs.hyb.sps;sps.hyb;hbs.bst},...
%         mm_,dlg(j_),...
%         {'';'';'PBS';'SP96';'HYB'},...
%         (cpp{j_}),bhr.nm{mm_},...
%         fnm,...
%         vtm_shift(2)*ones(5,1),...
%         {'none','none','none','none','none'},...
%         {':',':','-','-','-'});
%     % psa
%     set(0,'defaultaxescolororder',clr.rec_pbs_emp_hyb);
%     syn2ann_plot_compare(...
%         [1 0 0 0 0],...
%         {rec.org;pbs.hyb.sps;sps.hyb;hbs.bst},...
%         mm_,dlg(j_),...
%         {'REC';'PBS';'SP96';'HYB'},...
%         (cpp{j_}),bhr.nm{mm_},...
%         fnm,...
%         vtm_shift(2)*ones(5,1),...
%         {'none','none','none','none','none'},...
%         {'-',':','-.','-'});
%     
%     %% *POST-PROCESS - RECORDS vs SPECTRAL MATCHED (EMP)*
%     fnm = fullfile(sp,sprintf('%s_rec_spm_EMP_%s',st,ann.mtd.scl{j_}));
%     set(0,'defaultaxescolororder',clr.rec_spm);
%     syn2ann_plot_compare(...
%         [0 0 1 1 1],...
%         {rec.org;},...
%         mm_,dlg(j_),...
%         {'REC';'SPM'},...
%         (cpp{j_}),bhr.nm{mm_},...
%         fnm,...
%         [0,vtm_shift(1)]);
%     %% *POST-PROCESS - RECORDS vs PBS vs SPECTRAL MATCHED (EMP)*
%     fnm = fullfile(sp,sprintf('%s_rec_pbs_spm_EMP_%s',st,ann.mtd.scl{j_}));
%     set(0,'defaultaxescolororder',clr.rec_pbs_spm);
%     syn2ann_plot_compare(...
%         [1 1 0 0 0],...
%         {rec.org;pbs.bst;spm.sps.(cpp{j_})},...
%         mm_,dlg(j_),...
%         {'REC';'PBS';'SPM'},...
%         (cpp{j_}),bhr.nm{mm_},...
%         fnm,[0,vtm_shift(1),vtm_shift(1)],...
%         {'none';'none';'none'},...
%         {'-',':','-'});
end
