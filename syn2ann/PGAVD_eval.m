%% *Seismic Signals Peak Values*
% _Editor: Filippo Gatti
% CentraleSupélec - Laboratoire MSSMat
% Politecnico di Milano - DICA
% Copyright 2015_
%% NOTES
% _PGAVD_eval_: function to evaluate peak ground
% acceleration/velocity/displacement and time steps
%% INPUT:  
% * _dtm (sampling time step)_
% * _tha (input acceleration)_
% * _thv (input velocity)_
% * _thd (input displacement)_
%% OUTPUT: t_pga (time of pga)
% * _pga   (peak ground acceleration)_
% * _t_pgv (time of pgv)_
% * _pgv   (peak ground velocity)_
% * _t_pgd (time of pgd)_
% * _pgd   (peak ground displacement)_
function [varargout] = PGAVD_eval(varargin)
    
    %% SET-UP
    dtm = varargin{1};
    tha = varargin{2}(:);
    %%
    % _time vector_
    vtm = dtm*(0:(numel(tha)-1))';
    
    %% PEAK GROUND ACCELERATION
    [pga,idx] = max(abs(tha));
    vtm_pga   = vtm(idx);
    pga       = tha(idx);
    
    %% PEAK GROUND VELOCITY
    if nargin>=3
        thv=varargin{3};
    else
        thv=cumtrapz(tha)*dtm;
    end
    
    [pgv,idx] = max(abs(thv));
    vtm_pgv   = vtm(idx);
    pgv       = thv(idx);
    
    %% PEAK GROUND DISPLACEMENT
    if nargin==4
        thd = varargin{4};
    else
        thd = cumtrapz(thv)*dtm;
    end
    [pgd,idx] = max(abs(thd));
    vtm_pgd   = vtm(idx);
    pgd       = thd(idx);
    
    %% OUTPUT
    varargout{1} = vtm_pga;
    varargout{2} = pga;
    varargout{3} = vtm_pgv;
    varargout{4} = pgv;
    varargout{5} = vtm_pgd;
    varargout{6} = pgd;
    return
end

