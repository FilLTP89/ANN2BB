%% *Resampling seismic records*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _seismo_rsmpl_: function resample seismic records
%% INPUT:
% * _dtm (original time-step)_
% * _tha (original record)_
% * _fac (resampling rate)_
% * _scl (scale factor)_
%% OUTPUT:
% * _dtm_rsmpl (resampled time-step)_
% * _tha_rsmpl (resampled record)_
% * _ntm_rsmpl (resampled time-points)_
% * _vtm_rsmpl (resampled time-vector)_
function [varargout] = seismo_rsmpl(varargin)
    
    %% SET-UP
    
    dtm = varargin{1};
    tha = varargin{2}(:);
    fac = varargin{3};
    scl = 1.;
    if nargin>3
        scl = varargin{4};
    end
    ntm = numel(tha);
    vtm = dtm*(0:ntm-1)';
    
    %%
    % _time vectors_
    dtm_rsmpl = dtm/fac;
    ntm_rsmpl = round(ntm*fac);
    vtm_rsmpl = dtm_rsmpl*(0:ntm_rsmpl-1)';
    
    %% RESAMPLING
    tha_rsmpl = scl*interp1(vtm,tha,vtm_rsmpl,'spline');
    
    %% OUTPUT
    varargout{1} = dtm_rsmpl;
    varargout{2} = tha_rsmpl;
    varargout{3} = ntm_rsmpl;
    varargout{4} = vtm_rsmpl;
    return
end
