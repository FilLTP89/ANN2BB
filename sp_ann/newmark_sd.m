function [varargout] = newmark_sd(varargin)
    %===============
    % Newmark method for direct integration
    % Editor: Chiara Smerzini/Filippo Gatti
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % newmark_sd: function to compute acceleration/velocity/displacement of
    % a SDOF, by exploiting Newmark's method.
    % INPUT: ag (ground acceleration)
    %        dtm (time step for numerical integration)
    %        Tn (vector of natural period)
    %        zeta (damping ratio)
    % OUTPUT: ymax (maximum relative displacement of single dof)
    % N.B.: interpolation of input acceleration TH to improve accuracy of numerical
    % integration (step = dtm)
    %===============
    %======================================================================
    % SET-UP
    %======================================================================
    ag  = varargin{1};   % ag = interp1(t1,gacc1,t);
    dtm = varargin{2};
    Tn  = varargin{3};
    zeta = varargin{4};
    ntm = numel(ag);
    nTn = numel(Tn);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % newmark coefficients
    %----------------------------------------------------------------------
    beta = 0.25;
    gamma = 0.5;
    %----------------------------------------------------------------------
    % initialization
    %----------------------------------------------------------------------
    % spectral displacement
    sd  = -ones(nTn,1);
    % spectral velocity
    sv  = -ones(nTn,1);
    % spectral acceleration
    sa  = -ones(nTn,1);
    %======================================================================
    %======================================================================
    % NEWMARK INTEGRATION
    %======================================================================
    for j_ = 1:nTn % natural periods
        %------------------------------------------------------------------
        % natural circular frequency of sdof system
        %------------------------------------------------------------------
        wn = 2*pi/Tn(j_);
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % initial conditions
        %------------------------------------------------------------------
        y   = zeros(ntm,1);
        yp  = zeros(ntm,1);
        ypp = zeros(ntm,1);
        y0  = 0.;
        yp0 = 0.;
        y(1)   = y0;
        yp(1)  = yp0;
        ypp(1) = -ag(1)-2*wn*zeta*yp0-wn^2*y0;
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % Integration coefficients
        %------------------------------------------------------------------
        keff = wn^2 + 1/(beta*dtm^2) + gamma*2*wn*zeta/(beta*dtm);
        a1 = 1/(beta*dtm^2)+gamma*2*wn*zeta/(beta*dtm);
        a2 = 1/(beta*dtm)+2*wn*zeta*(gamma/beta-1);
        a3 = (1/(2*beta)-1)+2*wn*zeta*dtm*(gamma/(2*beta)-1);
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        % Newmark time scheme
        %------------------------------------------------------------------
        for i_ = 1:ntm-1 % time steps
            y(i_+1)   = (-ag(i_+1)+a1*y(i_)+a2*yp(i_)+a3*ypp(i_))/keff;
            ypp(i_+1) = ypp(i_)+...
                (y(i_+1)-y(i_)-dtm*yp(i_)-dtm^2*ypp(i_)/2)/(beta*dtm^2);
            yp(i_+1)  = yp(i_)+dtm*ypp(i_)+dtm*gamma*(ypp(i_+1)-ypp(i_));
        end
        %------------------------------------------------------------------
        sd(j_) = max(abs(y));
        sv(j_) = max(abs(yp));
        sa(j_) = max(abs(ypp));
    end
    % pseudo-spectral displacement
    psa = sd.*((2*pi./Tn).^2);
    % pseudo-spectral velocity
    psv = sd.*(2*pi./Tn);    
    %======================================================================
    varargout{1} = sd;
    varargout{2} = sv;
    varargout{3} = sa;
    varargout{4} = psa;
    varargout{5} = psv;
    
    return
end