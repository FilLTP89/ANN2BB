%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_plot_res_station_: function to plot ANN test results per station
%% *N.B.*
% Need for:
% _syn2ann_plot_compare.m_

global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
fprintf('%s-with HYBRID\n',bhr.nm{mm_});

%% *POST-PROCESS - RECORDS vs ANN*
for j_ = 1:numel(cpp.rec)
    flag.rec = seismo_dir_conversion(cpp.rec{j_});
    if any(strcmpi(flag.rec,flag.ann))
        fn = fullfile(sp,sprintf('%s_rec_ann_%s',st,tst.mtd.scl{kk_}));
        set(0,'defaultaxescolororder',clr21);
        syn2ann_plot_compare(...
            [1 0 0 0 0],...
            {rec;ann.tst{kk_}},...
            mm_,dlg(j_),...
            {'REC';'ANN'},...
            (cpp.rec{j_}),bhr.nm{mm_},...
            fn,...
            zeros(4,1),...
            {'none';'o'},...
            {'-';'none'});
    end
end