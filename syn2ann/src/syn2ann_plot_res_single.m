%% *GENERATION OF STRONG GROUND MOTION SIGNALS BY COUPLING PHYSICS-BASED ANALYSIS WITH ARTIFICIAL NEURAL NETWORKS*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_plot_res_single_: function to plot syn2ann results at each
% station for a single site class
%% *N.B.*
% Need for:
% _get_axis_tick.m,syn2ann_plot_res_station_sp96.m,
% syn2ann_plot_res_station_exsim.m,syn2ann_plot_res_station_both.m,
% syn2ann_fancy_plot.m_
fprintf('---------------------\n7. PLOTTING RESULTS\n---------------------\n');
%% *SET UP*
global pfg xlm xlb xtk ylm ylb ytk grd scl mrk tit utd max_val
cd(wd);close all;
if ~mon.map.flg
    syn2ann_setup_axes_common;
    
    for mm_ = 1:bhr.ns
        %% *SET-UP*
        syn2ann_setup_axes_single;
        switch lower(hybrid_type)
            case 'sp96'
                %% *SABETTA & PUGLIESE 1996*
                syn2ann_plot_res_station_sp96;
            case 'exsim'
                %% *EXSIM*
                syn2ann_plot_res_station_exsim;
        end
        %     syn2ann_fancy_plot;
    end
end