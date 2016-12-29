%% *Create ITACA record filename*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _itaca_monitor_name_: function to create itaca filenames
%% INPUT:
% * st (station id)
% * ev (recorded event)
% * dv (borehole device)
% * cp (motion direction)
% * rc (motion component)
% * ni (network)
%% OUTPUT: 
% * str(monitor file name)
function [varargout] = itaca_monitor_name(varargin)
    %% *SET-UP*
    st = varargin{1};
    ev = varargin{2};
    dv = varargin{3};
    cp = varargin{4};
    rc = varargin{5};
    ni = varargin{6};
    
    switch upper(rc)
        case 'A'
            rc = 'ACC';
        case 'V'
            rc = 'VEL';
        case 'D'
            rc = 'DIS';
    end
    %% *FILE NAME*
    str = strcat(upper(ni{1}),'.',upper(st),upper(dv),'..',upper(ni{2}),upper(cp),'.D.',...
        upper(ev),'.C.',upper(rc),'.ASC');
   
    if nargin>6
        pt  = varargin{7};  % absolute path
        str = fullfile(pt,str);
    end
    
    varargout{1} = str;
    return
end