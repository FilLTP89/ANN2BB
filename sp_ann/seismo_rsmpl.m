function [varargout] = seismo_rsmpl(varargin)
    %===============
    % Resampling seismic records
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % seismo_rsmpl: function resample seismic records
    % INPUT:  dtm (original time-step)
    %         tha (original record)
    %         fac (resampling rate)
    %         scl (scale factor)
    % OUTPUT: tha_rsmpl (resampled record)     
    %         dtm_rsmpl (resampled time-step)
    %         ntm_rsmpl (resampled time-points)
    %         vtm_rsmpl (resampled time-vector)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    dtm = varargin{1};
    tha = varargin{2}(:);
    fac = varargin{3};
    scl = 1.;
    if nargin>3
        scl = varargin{4};
    end
    ntm = numel(tha);
    vtm = dtm*(0:ntm-1)';
    %======================================================================
    %----------------------------------------------------------------------
    % time vectors
    %----------------------------------------------------------------------
    dtm_rsmpl = dt/fac;
    ntm_rsmpl = ntm*fac;
    vtm_rsmpl = dtm_rsmpl*(0:ntm_rsmpl-1)';
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % resampling
    %----------------------------------------------------------------------
    tha_rsmpl = scl*interp1(vtm,tha,vtm_rsmpl,'spline');
    %======================================================================
    varargout{1} = dtm_rsmpl;
    varargout{2} = tha_rsmpl;
    varargout{3} = ntm_rsmpl;
    varargout{4} = vtm_rsmpl;
    return
end
