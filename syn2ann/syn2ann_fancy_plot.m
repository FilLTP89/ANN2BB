%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_fancy_plot_: function to plot syn2ann results in a fancy way
%% *N.B.*
% Need for:
% _syn2ann_plot_aside.m_

global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
fprintf('%s-with HYBRID\n',bhr.nm{mm_});
dlg = cell(3,1);
dlg{1} = 'EW';
dlg{2} = 'NS';
dlg{3} = 'UD';
set(0,'defaultaxescolororder',clr0f);

for j_ = 1:numel(cpp)
    cpn.nm = cpp{j_};
    cpn.nb = j_;
    %% *POST-PROCESS - RECORDS VS NUMERICAL SIMULATIONS VS SPECTRAL-MATCHED*
    fnm = fullfile(sp,sprintf('%s_rec_sim_spm_%s',st,ann.mtd.scl{j_}));
    syn2ann_plot_aside(...
        ones(2,1),...
        {rec.org,nss.org,spm.sps.(cpp{j_})},...
        mm_,dlg(j_),...
        {'REC';'PBS';'SPM'},...
        cpn,bhr.nm{mm_},...
        fnm,...
        [0,vtm_shift(2),vtm_shift(2)]);
end