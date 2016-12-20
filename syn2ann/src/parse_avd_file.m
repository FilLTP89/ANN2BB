%% *Parse AVD record files*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _parse_avd_file: function to parse AVD formatted files
%% INPUT:
% * fn (file name)
%% OUTPUT:
% * _dtm (time-step)_
% * _ntm (number of time-steps)_
% * _vtm (time-vector)_
% * _tha (acceleration time-history)_
% * _thv (velocity time-history)_
% * _thd (displacement time-history)_

function [varargout] = parse_avd_file(varargin)
    %% *SET-UP*
    fn  = varargin{1};
    scl = varargin{2};
    
    data = importdata(fn);
    
    vtm = data(:,1);
    dtm = mean(diff(vtm));
    ntm = numel(vtm);
    tha = data(:,2)*scl;
    thv = data(:,3)*scl;
    thd = data(:,4)*scl;
    
    %% *OUTPUT*
    varargout{1} = dtm;
    varargout{2} = ntm;
    varargout{3} = vtm;
    varargout{4} = tha;
    varargout{5} = thv;
    varargout{6} = thd;
    return
end