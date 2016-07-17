% Sabetta-Pugliese synthetics accelerograms
% _Editor: Chiara Smerzini/Filippo Gatti
% DICA - Politecnico di Milano
% Copyright 2016_
%% NOTES
% _sabetta_: function for the generation of synthetics non-stationary accelerograms
% adapted from original fortran program of Sabetta F. and Pugliese A.
%% INPUT:
% * _mw (moment magnitude)_
% * _dep (epicentral distance [km])_
% * _scc (site conditions (0=rock; 1=shallow all.; 2=deep
%         alluvium))_
% * _sst (st. dev. of GMPE (0=median value,1=84th percentile)_
% * _dtm (time step of output accelerogram)_
% * _na  (number of accelerograms to be generated)_
% * _scl (scale factor [1=cm/s/s])_
% * _ivd (flag for output in velocity and displacement)_
%% OUTPUT:
% * _vtm (time vector)_
% * _tha (acceleration time histories [in cm/s/s])_
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
function [varargout] = sabetta(varargin)
    %% SET-UP
    if nargin > 1
        mw  = varargin{1};
        dep = varargin{2};
        scc = varargin{3}; % site conditions (0=rock; 1=shallow all.; 2=deep alluvium); 
        sst = varargin{4}; % st. dev. of GMPE (0=median value,1=84th percentile)
        dtm = varargin{5};
        scl = varargin{6};
    else
        mw  = varargin{1}.mw;
        dep = varargin{1}.dep;
        scc = varargin{1}.scc;
        sst = varargin{1}.sst;
        dtm = varargin{1}.dtm;
        scl = varargin{1}.scl;
    end
    %
    % _site conditions_
    %
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
    
    %% STRONG GROUND MOTION INDICATORS - EMPIRICAL ESTIMATION
    %
    % _strong ground motion duration
    % [Sabetta,Pugliese-1996],[Vanmarcke,Lai-1980]_
    %
    DV = 10^(-0.783 + 0.193*mw + 0.208*log10((dep^2 + 5.1^2)^(0.5)) -...
        0.133*S1 + 0.138*S2 + 0.247*sst);
    %
    % _arias intensity [Sabetta,Pugliese-1996]_
    %
    Ia = 10^( 0.729 + 0.911*mw - 1.818*log10((dep^2 + 5.3^2)^(0.5)) +...
        0.244*S1 + 0.139*S2 + 0.397*sst);
    %
    % _time delay [s] between S ans P waves (VP*VS/(VP-VS) = 7 km/s)_
    %
    T1 = dep/7;
    %
    % _other coefficients_
    %
    T2 = T1 + 0.5*DV;
    T3 = T1 + 2.5*DV; % = T2 + 2*DV
    TFc = T2 - 3.5 - dep/50;
    T_cost = T2 + 0.5*DV;
    %
    % _total duration of accelerogram_
    %
    tot_dur = 1.3*T3;
    
    T4 = tot_dur - T1;
    T_fond = T4/3;
    fo = 1/T_fond;
    
    %% INSTATANEOUS AVERAGE POWER
    %
    % _time vector_
    %
    vtm   = 0:dtm:tot_dur;
    ntm   = numel(vtm);
    t_val = zeros(ntm,1);
    t_val(1:ntm) = vtm - TFc;
    for i = 1:ntm
        if (t_val(i) < 1)
            t_val(i) = 1;
        end
        if (vtm(i) > T_cost)
            t_val(i) = t_val(i-1);
        end
    end
    %% FREQUENCY CONTENT
    %
    % _Nyquist frequency_
    fNy = 1/(2*dtm);
    %
    % _frequency vector_
    %
    vfr   = (fo:fo:fNy)';
    nfr   = numel(vfr);
    ind_f = (1:nfr)';
    %
    % _statistics - (NB: sqm_Pa.=2.5 in [Sabetta,Pugliese-1996] or =3 in
    % .for!)_
    sqm_Pa = log(T3/T2)/3;
    med_Pa = log(T2) + sqm_Pa^2;
    %
    % _Pa(t)_
    %
    Pa = Ia*lognpdf(vtm,med_Pa,sqm_Pa);
    %
    % _empirical regression for Fc [Hz] = central frequency_
    %
    Fc = exp(3.4 - 0.35.*log(t_val(:)) - 0.218*mw - 0.15*S2);
    %
    % _empirical regression for the ratio Fb/Fc ratio (frequency
    % bandwidth)_
    %
    Fb_Fc = 0.44 + 0.07*mw - 0.08*S1 + 0.03*S2;
    %
    % _statistics_
    %
    delta   = sqrt(log(1+Fb_Fc^2));
    ln_beta = log(Fc(:)) - 0.5*delta^2; 
    %% SYNTHETIC ACCELEROGRAMS
    tha = zeros(ntm,1);
    R = random('unif',0,2*pi,1,nfr)';
    Ccos = -ones(nfr,1);
    PS   = -ones(nfr,1);
    for i_ = 1:ntm
        %
        % PS in cm*cm/s/s/s
        %
        PS(:)   = (Pa(i_)./(ind_f(:).*sqrt(2*pi).*delta)).*...
            exp(-(log(vfr(:)) - ln_beta(i_)).^2./(2*delta^2));
        %
        % Ccos in cm/s/s
        %
        Ccos(:) = sqrt(2.*PS(:)).*cos(2.*pi.*vfr(:).*vtm(i_) + R(:));
        %
        % acc in cm/s/s
        %
        tha(i_) = sum(Ccos(:));
    end
    %
    % _scaling_
    %
    tha = tha.*scl;
    varargout{1} = vtm(:);
    varargout{2} = tha(:);
    return
end
