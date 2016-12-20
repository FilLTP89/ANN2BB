%% *Create KKNPP record filename*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _kknpp_monitor_name_: function to create kknpp filenames
%% INPUT:
% * st (station id)
% * ev (recorded event)
% * dv (borehole device)
% * cp (motion direction)
% * rc (motion component)
%% OUTPUT: 
% * str(monitor file name)
function [varargout] = kknpp_monitor_name(varargin)
    %% *SET-UP*
    st = varargin{1};
    ev = varargin{2};
    dv = varargin{3};
    %% *FILE NAME*
    str = strcat(upper(st),upper(dv),upper(ev),'.mat');
    
    if nargin>3
        pt  = varargin{4};  % absolute path
        str = fullfile(pt,str);
    end
    
    varargout{1} = str;
    return
end