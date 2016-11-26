%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_plot_res_compare_site_class_: function to plot tested ANN by comparing 
% site classes
%% *N.B.*
% Need for:
% _syn2ann_plot_compare.m_

global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd

st = bhr.nm{mm_};
ev = bhr.st{mm_}.ev{1};
fprintf('--------------------------\n');
fprintf('%s-%s\n',st,ev);

%% *POST-PROCESS - RECORDS vs ANN*
for j_ = 1:numel(cpp.rec)
    flag.rec = seismo_dir_conversion(cpp.rec{j_});
    if any(strcmpi(flag.rec,flag.ann))
        fnm = fullfile(spp,sprintf('%s_%s_rec_ann_cmp',st,ev));
        set(0,'defaultaxescolororder',clr21);
        
        syn2ann_plot_compare(...
            [1 0 0 0 0],...
            [{rec};ann.scp],...
            mm_,dlg(j_),...
            {'REC';'$T^\star-0.5 s$';'$T^\star-0.75 s$';'$T^\star-1.0 s$'},...
            (cpp.rec{j_}),bhr.nm{mm_},...
            fnm,...
            zeros(4,1),...
            {'none';'o';'d';'s'},...
            {'-';'--';':';'-.'});
    end
end
