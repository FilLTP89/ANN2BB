function [varargout] = speed_monitor_name(varargin)
    %===============
    % Create SPEED monitor filename
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % speed_monitor_name: function to create SPEED monitor filename from
    % working directory and monitor ID
    % INPUT:  id (monitor identity number)
    %         rc (motion component)
    % OUTPUT: str(monitor file name)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    id = varargin{1};
    rc = varargin{2};
    %======================================================================
    %======================================================================
    % FILE NAME COMPOSITION
    %======================================================================
    %----------------------------------------------------------------------
    % extra zeros
    %----------------------------------------------------------------------
    if id < 10
        ez ='0000';
    elseif id < 100   && id >= 10
        ez = '000';
    elseif id < 1000  && id >= 100
        ez = '00';
    elseif id < 10000 && id >= 1000
        ez = '0';
    else 
        ez = '';
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % file name
    %----------------------------------------------------------------------
    id = num2str(id);
    str = sprintf('monitor%s%s.%s',ez,id,rc);
    if nargin>2
        pt  = varargin{3};  % absolute path
        str = fullfile(pt,str);
    end
    %======================================================================
    varargout{1} = str;
    return
end