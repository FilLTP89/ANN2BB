%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _trann_test_plot_compare_: function to plot ANN test results at each
% station
%% *N.B.*
% Need for:
% _trann_plot_res_compare_siteclass.m_
fprintf('---------------------\n6. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd
cd(wd);close all;
% _color sets_
% records-pbs
clr0 = [rgb('Navy');rgb('IntenseBlue')];
% fil-fil-pbs-empirical/stochastic-hybrid
clr10 = [rgb('IntenseBlue');rgb('IntenseGreen');rgb('IntenseOrange')];
clr11 = [rgb('DarkGrey');rgb('DarkGrey');rgb('IntenseBlue');rgb('IntenseGreen');rgb('IntenseOrange')];
% records-pbs-empirical/stochastic-hybrid
clr12 = [rgb('Navy');rgb('IntenseBlue');rgb('IntenseGreen');rgb('IntenseOrange')];
% records-pbs-hybrid-hybrid
clr121 = [rgb('Orange');rgb('Black');[0.4 0.4 0.4];rgb('Black')];
% [rgb('Navy');rgb('IntenseBlue');rgb('Red');rgb('IntenseGreen')];
%[rgb('DarkIntenseOrange');rgb('Black');[0.5 0.5 0.5];[0.7 0.7 0.7]];
% record-hybrid-ann-spectral-matched
clr2 = [rgb('Navy');rgb('IntenseOrange');rgb('Red');rgb('Red')];
% records-spectral-matched
clr21 = [rgb('Navy');rgb('IntenseBlue');rgb('Red');rgb('IntenseOrange')];
% records-spectral-matched
clr3 = [rgb('Navy');rgb('Red')];
%cpp.rec = {'ew';'ns';'ud'};
%dlg = cell(3,1);
%dlg{1} = 'EW';
%dlg{2} = 'NS';
%dlg{3} = 'UD';

spp = '/home/filippo/Scrivania/ann_new';
if exist(spp,'dir')~=7
    spp = fullfile(filesep,'tmp1','gattif','heavy_images_new');
end
trann_setup_axes_common;
if size(cpp.ann)==1
    
    flag.ann = seismo_dir_conversion(cpp.ann);
    ann.scp = cell(tst.mtd.nr,1);
    for kk_=1:tst.mtd.nr
        ann.scp(kk_,1) = ann.tst(kk_);
    end
    for mm_ = 1:bhr.ns
        trann_setup_axes_single;
        disp('COMPARE BY NATURAL PERIOD!');
        trann_plot_res_compare_periods;
%         disp('COMPARE BY SITE CLASS!');
%         trann_plot_res_compare_siteclass;
    end
else
    disp('ERROR: ANN TRAINED ON DIFFERENT COMPONENTS!');
end
