%ANN2BB Broad-band strong ground motion generator coupling PBS and ANN
%
% ANN2BB_drive(flag_map,flag_plot_results)
%
% Multi-line paragraphs of descriptive text go here. It's fine for them to
% span lines. It's treated as preformatted text; help() and doc() will not
% re-wrap lines. In the editor, you can highlight paragraphs, right-click,
% and choose "Wrap selected comments" to re-flow the text.
%
% More detailed help is in the <a href="matlab: help foo>extended_help">extended help</a>.
% It's broken out like this so you can keep the main "help foo" text on
% a single screen, and then break out obscure parts to separate sections.
%
% Examples:
% foo(1,2,3)
%
% See also:
% BAR
% SOMECLASS/SOMEMETHOD

% *
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% *NOTES*
% _syn2ann_drive_: function to match the low-frequency synthetics' spectra
% from numerical simulations (SPEED/HISADA/SEM3D) to target spectra obtained via
% Artificial Neural Networks.
%% *N.B.*
% By setting flag_map=1, the user can run a fast analysis to plot single station
% time-histories/spectra (flag_map=0, demonstrative case) or run an analysis to
% obtain shake maps (flag_map=1, parallel implementation on remote cluster).
% Need for:
% _syn2ann_setup_maps.m,syn2ann_run_maps.m,
% syn2ann_setup_fast.m,syn2ann_run_fast.m_

function ANN2BB_drive(varargin)
    
    flag_map = 0; % flag to produce map data
    flag_plot_results =1; % flag to plot results
    if nargin>1
        flag_map=varargin{1};
    end
    if nargin==2
        flag_plot_results=varargin{2};
    end
    
    if flag_map % write map data
        %% *1). CUSTOMIZE ANALYSIS SET-UP*
        syn2ann_setup_maps;
        
        %% *2). RUN SYN2ANN TO GET SHAKE MAPS (DNC)*
        syn2ann_run_maps;
    else
        
        %% *1). CUSTOMIZE ANALYSIS SET-UP*
        syn2ann_setup_fast;
        
        %% *2). PARSING REC (DNC)*
        syn2ann_rec_drive;
        
        %% *3). RUN SYN2ANN ON SINGLE STATIONS (DNC)*
        syn2ann_run;
        
        %% *5). SAVE RESULTS (DNC)*
        syn2ann_save_res;
        
        %% *6). PLOT RESULTS*
        if flag_plot_results
            syn2ann_plot_res_single;
        end
    end
    return
end