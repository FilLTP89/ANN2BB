%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_plot_res_compare_periods_: function to plot tested ANN by comparing
% different training periods
%% *N.B.*
% Need for:
% _syn2ann_plot_compare.m_

global pfg xlm xlb xtk ylm ylb ytk grd scl mrk utd tit


%% *POST-PROCESS - RECORDS vs ANN*
for j_ = 1:numel(cpp.rec)
    flag.rec = seismo_dir_conversion(cpp.rec{j_});
    if any(strcmpi(flag.rec,flag.ann))
        fnm = fullfile(spp,sprintf('%s_%s_rec_ann_%s_%s_cpc',...
            st,fgn,ann.scp{1}.cpp,ann.scp{1}.scl));
        set(0,'defaultaxescolororder',clr121);
        
        %         syn2ann_plot_compare(...
        %             [1 0 0 0 0],...
        %             [{rec};ann.scp],...
        %             mm_,dlg(j_),...
        %             {'REC';'$T^\star-0.5 s$';'$T^\star-0.75 s$';'$T^\star-1.0 s$'},...
        %             (cpp.rec{j_}),ttt,...
        %             fnm,...
        %             zeros(4,1),...
        %             {'none';'none';'s';'v'},...
        %             {'-';'-';'-';'--'},...
        %             [6,2,0.8,1.5]);
        
    end
end
trann_test_psa_performance(rec.org,ann.scp,mm_,ttt,fnm);