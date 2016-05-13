function [varargout] = sabetta(varargin)
    %===============
    % Sabetta-Pugliese synthetics accelerograms
    % Editor: Chiara Smerzini/Filippo Gatti
    % DICA - Politecnico di Milano
    % Copyright 2016
    % NOTES
    % function for the generation of synthetics non-stationary accelerograms
    % adapted from original fortran program of Sabetta F. and Pugliese A.
    % INPUT:  mw (moment magnitude)
    %         dep (epicentral distance [km])
    %         scc (site conditions (0=rock; 1=shallow all.; 2=deep alluvium))
    %         sst (st. dev. of GMPE (0=median value,1=84th percentile)
    %         dtm (time step of output accelerogram)
    %         na  (number of accelerograms to be generated)
    %         scl (scale factor [1=cm/s/s])
    %         ivd (flag for output in velocity and displacement)
    % OUTPUT: vtm (time vector)
    %         acc (acceleration time histories [in cm/s/s])
    % REFERENCES:
    % @Article{Sabetta_Pugliese_1996,
    %   author =  {Sabetta, F. and Pugliese, A.},
    %   title =   {{Estimation of Response Spectra and Simulation of Nonstationary Earthquake Ground Motions}},
    %   journal = {Bulletin of the Seismological Society of America},
    %   year =    {1996},
    %   volume =  {86},
    %   number =  {2},
    %   pages =   {337--352}
    % }
    %===============
    
    %======================================================================
    % SET-UP
    %======================================================================
    if nargin > 1
        mw  = varargin{1};
        dep = varargin{2};
        scc = varargin{3};
        sst = varargin{4};
        dtm = varargin{5};
        na  = varargin{6};
        scl = varargin{7};
    else
        mw  = varargin{1}.mw;
        dep = varargin{1}.dep;
        scc = varargin{1}.scc;
        sst = varargin{1}.sst;
        dtm = varargin{1}.dtm;
        na  = varargin{1}.na;
        scl = varargin{1}.scl;
    end
    %======================================================================
    %======================================================================
    % SITE CONDITIONS
    %======================================================================
    switch scc
        case 0
            S1=0;
            S2=0;
        case 1
            S1=1;
            S2=0;
        case 2
            S1=0;
            S2=1;
    end
    %======================================================================
    %======================================================================
    % STRONG GROUND MOTION INDICATORS - EMPIRICAL ESTIMATION
    %======================================================================
    %----------------------------------------------------------------------
    % strong ground motion duration [Sabetta,Pugliese-1996],[Vanmarcke,Lai-1980]
    %----------------------------------------------------------------------
    DV = 10^(-0.783 + 0.193*mw + 0.208*log10((dep^2 + 5.1^2)^(0.5)) -...
        0.133*S1 + 0.138*S2 + 0.247*sst);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % arias intensity [Sabetta,Pugliese-1996]
    %----------------------------------------------------------------------
    Ia = 10^( 0.729 + 0.911*mw - 1.818*log10((dep^2 + 5.3^2)^(0.5)) +...
        0.244*S1 + 0.139*S2 + 0.397*sst);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % time delay [s] between S ans P waves (VP*VS/(VP-VS) = 7 km/s)
    %----------------------------------------------------------------------
    T1 = dep/7;
    %----------------------------------------------------------------------
    % other coefficients
    %----------------------------------------------------------------------
    T2 = T1 + 0.5*DV;
    T3 = T1 + 2.5*DV; % = T2 + 2*DV
    TFc = T2 - 3.5 - dep/50;
    T_cost = T2 + 0.5*DV;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % total duration of accelerogram
    %----------------------------------------------------------------------
    tot_dur = 1.3*T3;
    %----------------------------------------------------------------------
    T4 = tot_dur - T1;
    T_fond = T4/3;
    fo = 1/T_fond;
    %======================================================================
    %======================================================================
    % Pa(t) INSTATANEOUS AVERAGE POWER
    %======================================================================
    % Pa [cm^2/s^4] lognormal time envelope function depicting the spectral
    % amplitude variation of ground motion in time.
    %----------------------------------------------------------------------
    % time vector
    %----------------------------------------------------------------------
    vtm   = 0:dtm:tot_dur;
    ntm   = numel(vtm);
    t_val = zeros(1,ntm);
    t_val(1:ntm) = vtm - TFc;
    for i = 1:ntm
        if (t_val(i) < 1)
            t_val(i) = 1;
        end
        if (vtm(i) > T_cost)
            t_val(i) = t_val(i-1);
        end
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Nyquist frequency
    %----------------------------------------------------------------------
    fmax = 1/(2*dtm);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % statistics - (NB: sqm_Pa.=2.5 in [Sabetta,Pugliese-1996] or =3 in .for!)
    %----------------------------------------------------------------------
    sqm_Pa = log(T3/T2)/3;
    med_Pa = log(T2) + sqm_Pa^2;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Pa(t)
    %----------------------------------------------------------------------
    Pa = Ia*lognpdf(vtm,med_Pa,sqm_Pa);
    %======================================================================
    %======================================================================
    % FREQUENCY CONTENT
    %======================================================================
    %----------------------------------------------------------------------
    % empirical regression for Fc [Hz] = central frequency
    %----------------------------------------------------------------------
    Fc = exp(3.4 - 0.35.*log(t_val) - 0.218*mw - 0.15*S2);
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % empirical regression for the ratio Fb/Fc ratio (frequency bandwidth)
    %----------------------------------------------------------------------
    Fb_Fc = 0.44 + 0.07*mw - 0.08*S1 + 0.03*S2;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % statistics
    %----------------------------------------------------------------------
    delta   = sqrt(log(1+Fb_Fc^2));
    ln_beta = log(Fc) - 0.5*delta^2;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % frequency vector
    %----------------------------------------------------------------------
    frq   = fo:fo:fmax;
    nfr   = length(frq);
    ind_f = 1:nfr;
    %======================================================================
    %======================================================================
    % SYNTHETIC ACCELEROGRAMS
    %======================================================================
    a = zeros(ntm,na);
    % Ccos_vel = zeros(1,nfr);
    % Ccos_dis = zeros(1,nfr);
    % vel = zeros(nt,na);
    % dis = zeros(nt,na);
    for k_ = 1:na % number of synthetics
        R=random('unif',0,2*pi,1,nfr);
        for i_=1:ntm
            % PS in cm^2 / s^3
            PS  = (Pa(i_)./(ind_f.*sqrt(2*pi).*delta)).*exp(-(log(frq) -...
                ln_beta(i_)).^2./(2*delta^2));
            % Ccos in cm / s^2
            Ccos = sqrt(2.*PS).*cos(2.*pi.*frq.*vtm(i_) + R);
            %         % Ccos_vel in cm / s
            %         Ccos_vel(1:nfr) = Ccos(1:nfr)./(2*pi.*f);
            %         % Ccos_dis in cm
            %          Ccos_dis(1:nfr) = Ccos(1:nfr)./(2*pi.*f).^2;
            % acc in cm/s/s
            a(i_,k_) = sum(Ccos);
            %         % vel in cm/s
            %         vel(i,k) = sum(Ccos_vel(1:nfr));
            %         % dis in cm
            %         dis(i,k) = sum(Ccos_dis(1:nfr));
        end
    end
    %----------------------------------------------------------------------
    % scaling
    %----------------------------------------------------------------------
    a = a.*scl;
    %======================================================================
    
    varargout{1} = vtm(:);
    varargout{2} = a;
    return
end
