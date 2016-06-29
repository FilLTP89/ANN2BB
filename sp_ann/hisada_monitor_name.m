function [varargout] = hisada_monitor_name(varargin)
    %===============
    % Create HISADA monitor filename
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % hisada_monitor_name: function to create HISADA monitor filename from
    % working directory and monitor ID
    % INPUT:  id (monitor identity number)
    %         cp (motion axis)
    % OUTPUT: str(monitor file name)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    id = varargin{1};
    cp = varargin{2};
    
    switch(lower(cp))
        case 'e'
            cp = 'y';
        case 'n'
            cp = 'x';
    end
    %======================================================================
    %======================================================================
    % FILE NAME COMPOSITION
    %======================================================================
    %----------------------------------------------------------------------
    % file name
    %----------------------------------------------------------------------
    str = sprintf('%s_dat.%u',cp,id);
    if nargin>2
        pt  = varargin{3};      % absolute path
        str = fullfile(pt,str);
    end
    %======================================================================
    varargout{1} = str;
    return
end