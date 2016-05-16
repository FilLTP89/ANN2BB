%% *Newmark method for direct integration*
% _Editor: Chiara Smerzini/Filippo Gatti
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _newmark_sd_: function to compute acceleration/velocity/displacement of
% a SDOF, by exploiting Newmark's method.
%% INPUT: 
% * _ag (ground acceleration)_
% * _dtm (time step for numerical integration)_
% * _vTn (vector of natural period)_
% * _zeta (damping ratio)_
%% OUTPUT: 
% * _ymax (maximum relative displacement of single dof)_
%% N.B.: 
% interpolation of input acceleration TH to improve accuracy of numerical
% integration (step = dtm)
function [varargout] = newmark_sd(varargin)
    
    %% SET-UP
    ag  = varargin{1};   % ag = interp1(t1,gacc1,t);
    dtm = varargin{2};
    vTn  = varargin{3};
    zeta = varargin{4};
    ntm = numel(ag);
    nTn = numel(vTn);
    
    %%
    % _newmark coefficients_
    beta = 0.25;
    gamma = 0.5;
    %%
    % _spectral displacement_
    sd  = -ones(nTn,1);
    %%
    % _spectral velocity_
    sv  = -ones(nTn,1);
    %%
    % _spectral acceleration_
    sa  = -ones(nTn,1);
    
    %% NEWMARK INTEGRATION
    
    for j_ = 1:nTn % natural periods
        %%
        % _natural circular frequency of sdof system_
        wn = 2*pi/vTn(j_);
        %%
        % _initial conditions_
        y   = zeros(ntm,1);
        yp  = zeros(ntm,1);
        ypp = zeros(ntm,1);
        y0  = 0.;
        yp0 = 0.;
        y(1)   = y0;
        yp(1)  = yp0;
        ypp(1) = -ag(1)-2*wn*zeta*yp0-wn^2*y0;
        %%
        % _integration coefficients_
        keff = wn^2 + 1/(beta*dtm^2) + gamma*2*wn*zeta/(beta*dtm);
        a1 = 1/(beta*dtm^2)+gamma*2*wn*zeta/(beta*dtm);
        a2 = 1/(beta*dtm)+2*wn*zeta*(gamma/beta-1);
        a3 = (1/(2*beta)-1)+2*wn*zeta*dtm*(gamma/(2*beta)-1);
        %%
        % _Newmark time scheme_
        for i_ = 1:ntm-1 % time steps
            y(i_+1)   = (-ag(i_+1)+a1*y(i_)+a2*yp(i_)+a3*ypp(i_))/keff;
            ypp(i_+1) = ypp(i_)+...
                (y(i_+1)-y(i_)-dtm*yp(i_)-dtm^2*ypp(i_)/2)/(beta*dtm^2);
            yp(i_+1)  = yp(i_)+dtm*ypp(i_)+dtm*gamma*(ypp(i_+1)-ypp(i_));
        end
        %%
        % _maximum values_
        sd(j_) = max(abs(y));
        sv(j_) = max(abs(yp));
        sa(j_) = max(abs(ypp));
    end
    %%
    % _pseudo-spectral displacement_
    psa = sd.*((2*pi./vTn).^2);
    %%
    % _pseudo-spectral velocity_
    psv = sd.*(2*pi./vTn);
    %% OUTPUT
    varargout{1} = sd;
    varargout{2} = sv;
    varargout{3} = sa;
    varargout{4} = psa;
    varargout{5} = psv;
    
    return
end