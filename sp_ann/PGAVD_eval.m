function [varargout] = PGAVD_eval(varargin)
    %===============
    % Seismic Signals' Peak Values
    % Editor: Filippo Gatti
    % CentraleSupélec - Laboratoire MSSMat
    % Politecnico di Milano - DICA
    % Copyright 2015
    % NOTES
    % PGAVD_eval: function to evaluate peak ground
    % acceleration/velocity/displacement and time steps
    % INPUT:  dtm (sampling time step)
    %         tha (input acceleration)
    %         thv (input velocity)
    %         thd (input displacement)
    % OUTPUT: t_pga (time of pga)
    %         pga   (peak ground acceleration)
    %         t_pgv (time of pgv)
    %         pgv   (peak ground velocity)
    %         t_pgd (time of pgd)
    %         pgd   (peak ground displacement)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    dtm = varargin{1};
    tha = varargin{2}(:);
    
    vtm = dtm*(0:(numel(tha)-1))';
    %======================================================================
    %======================================================================
    % PEAK GROUND ACCELERATION
    %======================================================================
    [pga,idx] = max(abs(tha));
    vtm_pga     = vtm(idx);
    
    varargout{1} = vtm_pga;
    varargout{2} = pga;
    %======================================================================
    %======================================================================
    % PEAK GROUND VELOCITY
    %======================================================================
    if nargin>=3
        thv=varargin{3};
    else
        thv=cumtrapz(tha)*dtm;
    end
    
    [pgv,idx] = max(abs(thv));
    vtm_pgv   = vtm(idx);
    
    varargout{3} = vtm_pgv;
    varargout{4} = pgv;
    %======================================================================
    %======================================================================
    % PEAK GROUND DISPLACEMENT
    %======================================================================
    if nargin==4
        thd = varargin{4};
    else
        thd = cumtrapz(thv)*dtm;
    end
    [pgd,idx] = max(abs(thd));
    vtm_pgd   = vtm(idx);
    %======================================================================
    varargout{1} = vtm_pga;
    varargout{2} = pga;
    varargout{3} = vtm_pgv;
    varargout{4} = pgv;
    varargout{5} = vtm_pgd;
    varargout{6} = pgd;
    return
end

